#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mdtraj as md
import argparse
import numpy as np
import sys
import os
import glob

def identify_pocket_atoms(ref_traj, cutoff=1.0):
    """Identify protein backbone pocket atoms within a given distance of the ligand.
    If ligand is not present, use all protein backbone heavy atoms."""
    # Select protein backbone (chain A) and ligand (not protein, no H)
    protein_atoms = ref_traj.topology.select("protein and backbone and not element H")
    ligand_atoms = ref_traj.topology.select("not protein and not element H")

    if len(protein_atoms) == 0:
        print("❌ Error: No protein backbone atoms found in the reference PDB.")
        sys.exit(1)

    # Fallback for protein-only case
    if len(ligand_atoms) == 0:
        print("⚠ Warning: No ligand found. Using protein backbone heavy atoms as pocket.")
        return protein_atoms, False

    # Compute distances between protein atoms and ligand atoms
    pairs = np.array([[p, l] for p in protein_atoms for l in ligand_atoms])
    distances = md.compute_distances(ref_traj, pairs).reshape(len(protein_atoms), len(ligand_atoms))

    # Identify pocket atoms (within cutoff distance)
    pocket_mask = np.any(distances < cutoff, axis=1)
    pocket_atoms = protein_atoms[pocket_mask]

    if len(pocket_atoms) == 0:
        print("⚠ Warning: No pocket residues found in reference. Using full backbone instead.")
        return protein_atoms, True

    return pocket_atoms, True

def get_matched_pocket_atoms(ref_traj, target_traj, ref_pocket_atoms):
    target_lookup = {}
    for atom in target_traj.topology.atoms:
        key = (atom.name, atom.residue.name, atom.residue.index, atom.residue.chain.index)
        if key not in target_lookup:
            target_lookup[key] = atom.index  # Save first match only

    matched_target_atoms = []
    for idx in ref_pocket_atoms:
        ref_atom = ref_traj.topology.atom(idx)
        key = (ref_atom.name, ref_atom.residue.name, ref_atom.residue.index, ref_atom.residue.chain.index)
        if key in target_lookup:
            matched_target_atoms.append(target_lookup[key])

    return matched_target_atoms

def main():
    parser = argparse.ArgumentParser(description="Align a PDB file to a reference using the protein backbone pocket residues (MDTraj).")
    parser.add_argument("target", help="Target PDB file to be aligned.")
    parser.add_argument("--ref", default="model.pdb", help="Reference PDB file (default: model.pdb)")
    parser.add_argument("--cutoff", type=float, default=1.0, help="Pocket cutoff distance in nm (default: 1.0)")
    parser.add_argument("-o", "--output", default="aligned_model.pdb", help="An aligned model filename (default: aligned_model.pdb at the same input directory)")
    parser.add_argument("--oligand", default="aligned_ligand.pdb", help="An aligned ligand filename (default: aligned_ligand.pdb at the same input directory)")
    args = parser.parse_args()

    print(f"📂 Loading reference: {args.ref}")
    pattern = os.path.expanduser(args.ref)
    ref_pdb_file = glob.glob(pattern)
    ref_pdb_file = os.path.abspath(ref_pdb_file[0])
    print(f"  - {ref_pdb_file}")
    ref_traj = md.load(args.ref)
    print(ref_traj)  # Print summary of the reference trajectory

    print(f"📂 Loading target: {args.target}")
    pattern = os.path.expanduser(args.target)
    target_pdb_files = glob.glob(pattern)
    target_pdb_files = [os.path.abspath(f) for f in target_pdb_files]
    target_pdb_files.sort(key=lambda x: (not os.path.basename(os.path.dirname(x)).startswith("best_pose"), x))
    print("\n".join(f"  - {f}" for f in target_pdb_files))
    target_traj = md.load(target_pdb_files)
    print(target_traj)  # Print summary of the target trajectory

    # Determine output directory (same as reference PDB?)
    if args.output == "aligned_model.pdb":
        output_files = [os.path.join(os.path.dirname(pdb_path), "aligned_model.pdb") for pdb_path in target_pdb_files]
    else:
        output_files = args.output

    # Determine oligand directory (same as reference PDB?)
    if args.oligand == "aligned_ligand.pdb":
        oligand_files = [os.path.join(os.path.dirname(pdb_path), "aligned_ligand.pdb") for pdb_path in target_pdb_files]
    else:
        oligand_files = args.oligand

    print("🔍 Identifying protein backbone pocket atoms from reference...")
    # Add flag to indicate if ligand exists
    ref_pocket_atoms, has_ligand = identify_pocket_atoms(ref_traj, cutoff=args.cutoff)
    target_pocket_atoms = get_matched_pocket_atoms(ref_traj, target_traj, ref_pocket_atoms)

    if len(ref_pocket_atoms) == len(target_pocket_atoms):
        print("📐 Aligning target to reference using protein backbone pocket atoms...")
        target_traj.superpose(ref_traj, atom_indices=target_pocket_atoms, ref_atom_indices=ref_pocket_atoms)

        # Compute RMSD of protein backbone
        idx = target_traj.topology.select("protein and backbone and not element H")
        target_protein_traj = target_traj.atom_slice(idx)
        idx = ref_traj.topology.select("protein and backbone and not element H")
        ref_protein_traj = ref_traj.atom_slice(idx)
        rmsd_values = md.rmsd(target_protein_traj, ref_protein_traj)
        print(f"📊 Pocket-Aligned Protein Backbone RMSD (MDTraj): {rmsd_values} nm")

        # Compute RMSD of protein pocket
        target_pocket_traj = target_traj.atom_slice(target_pocket_atoms)
        ref_pocket_traj = ref_traj.atom_slice(ref_pocket_atoms)
        rmsd_values = md.rmsd(target_pocket_traj, ref_pocket_traj)
        print(f"📊 Pocket-Aligned Protein Backbone Pocket RMSD (MDTraj): {rmsd_values} nm")

        # Save aligned model structures
        print(f"💾 Saving aligned model structures")
        if args.output == "aligned_model.pdb":
            for i in range(len(target_traj)):
                print(f" - {output_files[i]}")
                target_traj[i].save(output_files[i])
        else:
            print(f" - {output_files}")
            target_traj.save(output_files)

        # Only run ligand saving if ligand is detected
        if has_ligand:
            idx = target_traj.topology.select("not protein and not element H")
            target_ligand_traj = target_traj.atom_slice(idx)
            print(f"💾 Saving aligned ligand structures")
            if args.oligand == "aligned_ligand.pdb":
                for i in range(len(target_ligand_traj)):
                    print(f" - {oligand_files[i]}")
                    target_ligand_traj[i].save(oligand_files[i])
            else:
                print(f" - {oligand_files}")
                target_ligand_traj.save(oligand_files)
        else:
            print("ℹ️ No ligand found — skipping ligand output.")

        print("✅ Alignment complete.\n")

if __name__ == "__main__":
    main()

