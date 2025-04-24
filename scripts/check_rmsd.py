import argparse
import mdtraj as md
import numpy as np
import glob
import os

def is_same_traj(ref_traj, target_traj):
    """
    Compare atom info between two MDTraj slices (atom name, residue name, index, chain ID).
    Returns (True, None) if matched, or (False, first_mismatch_info) if not.
    """
    ref_atoms = ref_traj.topology.atoms
    tgt_atoms = target_traj.topology.atoms

    if target_traj.n_atoms != ref_traj.n_atoms:
        return False, ("length mismatch", target_traj.n_atoms, ref_traj.n_atoms)

    for i, (a1, a2) in enumerate(zip(ref_atoms, tgt_atoms)):
        ref_info = (a1.name, a1.residue.name, a1.residue.index, a1.residue.chain.index)
        tgt_info = (a2.name, a2.residue.name, a2.residue.index, a2.residue.chain.index)
        if ref_info != tgt_info:
            return False, (i, ref_info, tgt_info)

    return True, None


def compute_rmsd(ref_slice, target_slice, description):
    """Compute RMSD using MDTraj and a direct method, ensuring atom names match."""

    if ref_slice.n_atoms == 0:
        print(f"❌ Error: No {description} atoms found in the reference PDB. Skipping RMSD calculation.")
        return

    if target_slice.n_atoms == 0:
        print(f"❌ Error: No {description} atoms found in the target PDB. Skipping RMSD calculation.")
        return

    # Check the same atom info before RMSD calculations
    same, mismatch = is_same_traj(ref_slice, target_slice)
    if same:
        # Direct RMSD
        diff = target_slice.xyz - ref_slice.xyz[0]
        rmsd_direct = np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=1))
        print(f"📊 {description} RMSD (Direct): Min = {rmsd_direct.min():.3f} nm ; {rmsd_direct} nm")

        # MDTraj RMSD
        rmsd_mdtraj = md.rmsd(target_slice, ref_slice)
        print(f"📊 {description} RMSD (MDTraj): Min = {rmsd_mdtraj.min():.3f} nm ; {rmsd_mdtraj} nm")
    else:
        print(f"⚠ Warning: Mismatch found between reference and target for {description}.")
        if mismatch[0] == "length mismatch":
            print(f"   Atom count differs: reference={mismatch[1]}, target={mismatch[2]}")
        else:
            i, ref_atom, tgt_atom = mismatch
            print(f"   First mismatch at the atom info (atom_name, residue_name, residue_index, chain_index):")
            print(f"     Reference: {ref_atom}")
            print(f"     Target:    {tgt_atom}")
        print(f"   Skipping RMSD calculation for {description}.")

def identify_pocket_atoms(ref_traj, cutoff=1.0):
    """Identify protein backbone pocket atoms within a given distance of the ligand."""
    protein_atoms = ref_traj.topology.select("protein and backbone and not element H")
    ligand_atoms = ref_traj.topology.select("not protein and not element H")

    if len(protein_atoms) == 0:
        return np.array([])

    if len(ligand_atoms) == 0:
        return np.array([])

    # Compute distances between protein and ligand atoms
    distances = md.compute_distances(ref_traj, np.array(np.meshgrid(protein_atoms, ligand_atoms)).T.reshape(-1, 2))
    distances = distances.reshape(len(protein_atoms), len(ligand_atoms))

    # Identify pocket residues within the cutoff distance
    pocket_atoms = protein_atoms[np.any(distances < cutoff, axis=1)]

    if len(pocket_atoms) == 0:
        print(f"⚠ Warning: No pocket residues found within {cutoff} nm.")

    return pocket_atoms

def get_matched_pocket_atoms(ref_traj, target_traj, ref_pocket_atoms):
    # Step 1: Build a lookup table for atoms in the target_traj
    target_lookup = {}
    for atom in target_traj.topology.atoms:
        key = (atom.name, atom.residue.name, atom.residue.index, atom.residue.chain.index)
        if key not in target_lookup:
            target_lookup[key] = atom.index  # Save first match only

    # Step 2: Match each ref pocket atom using the lookup
    matched_target_atoms = []
    for idx in ref_pocket_atoms:
        ref_atom = ref_traj.topology.atom(idx)
        key = (ref_atom.name, ref_atom.residue.name, ref_atom.residue.index, ref_atom.residue.chain.index)

        if key in target_lookup:
            matched_target_atoms.append(target_lookup[key])

    return matched_target_atoms


def main():
    parser = argparse.ArgumentParser(description="Check RMSD between two PDB structures.")
    parser.add_argument("target_pdb", type=str, help="Path to the target PDB file")
    parser.add_argument("--ref", type=str, help="Path to the reference PDB file")
    parser.add_argument("--cutoff", type=float, default=1.0, help="Cutoff distance for pocket residue identification (default: 1.0 nm)")
    args = parser.parse_args()

    print(f"📂 Loading reference PDB: {args.ref}")
    try:
        pattern = os.path.expanduser(args.ref)
        ref_pdb_file = glob.glob(pattern)
        ref_pdb_file = os.path.abspath(ref_pdb_file[0])
        print(f"  - {ref_pdb_file}")
        ref_traj = md.load(ref_pdb_file)
        if ref_traj.n_frames == 0:
            raise ValueError("Reference PDB is empty.")
    except Exception as e:
        print(f"❌ Error: Failed to load reference PDB ({args.ref}): {e}")
        return

    print(f"📂 Loading target PDB: {args.target_pdb}")
    try:
        pattern = os.path.expanduser(args.target_pdb)
        pdb_files = glob.glob(pattern)
        pdb_files = [os.path.abspath(f) for f in pdb_files]
        pdb_files.sort(key=lambda x: (not os.path.basename(os.path.dirname(x)).startswith("best_pose"), x))
        print("\n".join(f"  - {f}" for f in pdb_files))
        target_traj = md.load(pdb_files)

        if target_traj.n_frames == 0:
            raise ValueError("Target PDB is empty.")
    except Exception as e:
        print(f"❌ Error: Failed to load target PDB ({args.target_pdb}): {e}")
        return

    # Checking the RMSD values
    description = "Protein Backbone"
    ref_backbone_idx = ref_traj.topology.select("protein and backbone and not element H")
    target_backbone_idx = target_traj.topology.select("protein and backbone and not element H")
    if len(ref_backbone_idx) > 0 and len(target_backbone_idx) > 0:
        ref_backbone = ref_traj.atom_slice(ref_backbone_idx)
        target_backbone = target_traj.atom_slice(target_backbone_idx)
        compute_rmsd(ref_backbone, target_backbone, description)
    else:
        print("⚠ Warning: Could not identify protein backbone in the target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")

    description = "Ligand"
    ref_ligand_idx = ref_traj.topology.select("not protein and not element H")
    target_ligand_idx = target_traj.topology.select("not protein and not element H")
    if len(ref_ligand_idx) > 0 and len(target_ligand_idx) > 0:
        ref_ligand = ref_traj.atom_slice(ref_ligand_idx)
        target_ligand = target_traj.atom_slice(target_ligand_idx)
        compute_rmsd(ref_ligand, target_ligand, description)
    else:
        print("⚠ Warning: Could not identify ligand in the target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")

    description = "Protein Backbone Pocket"
    ref_pocket_atoms = identify_pocket_atoms(ref_traj, args.cutoff)
    target_pocket_atoms = get_matched_pocket_atoms(ref_traj, target_traj, ref_pocket_atoms)
    if len(ref_pocket_atoms) > 0 and len(target_pocket_atoms) > 0:
        ref_pocket = ref_traj.atom_slice(ref_pocket_atoms)
        target_pocket = target_traj.atom_slice(target_pocket_atoms) # comparing the same pocket
        compute_rmsd(ref_pocket, target_pocket, description)
    else:
        print("⚠ Warning: No pocket residues found in target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")


if __name__ == "__main__":
    main()
