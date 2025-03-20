import argparse
import mdtraj as md
import numpy as np

def compute_rmsd(ref_slice, target_slice, description):
    """Compute RMSD using MDTraj and a direct method, ensuring atom names match."""

    if ref_slice.n_atoms == 0:
        print(f"❌ Error: No {description} atoms found in the reference PDB. Skipping RMSD calculation.")
        return

    if target_slice.n_atoms == 0:
        print(f"❌ Error: No {description} atoms found in the target PDB. Skipping RMSD calculation.")
        return

    # Extract atom names from topology
    ref_atom_names = [atom.name for atom in ref_slice.topology.atoms]
    target_atom_names = [atom.name for atom in target_slice.topology.atoms]

    # Ensure atom names are identical in both structures
    if ref_atom_names != target_atom_names:
        print(f"⚠ Warning: Atom names do not match between reference and target for {description}.")
        print(f"   Reference atoms: {ref_atom_names[:5]} ... ({len(ref_atom_names)} total)")
        print(f"   Target atoms:    {target_atom_names[:5]} ... ({len(target_atom_names)} total)")
        print(f"   Skipping RMSD calculation for {description}.")
        return

    # Direct RMSD calculation
    diff = target_slice.xyz - ref_slice.xyz[0]
    rmsd_direct = np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=1))
    print(f"📊 {description} RMSD (Direct): Min = {rmsd_direct.min():.3f} nm ; {rmsd_direct} nm")

    # MDTraj RMSD
    rmsd_mdtraj = md.rmsd(target_slice, ref_slice)
    print(f"📊 {description} RMSD (MDTraj): Min = {rmsd_mdtraj.min():.3f} nm ; {rmsd_mdtraj} nm")

def identify_pocket_residues(ref_traj, cutoff=1.0):
    """Identify protein backbone residues within a given distance of the ligand."""
    protein_atoms = ref_traj.topology.select("protein and backbone and not element H")
    ligand_atoms = ref_traj.topology.select("not protein and not element H")

    if len(protein_atoms) == 0:
        print(f"❌ Error: Could not identify protein in the reference structure.")
        return np.array([])

    if len(ligand_atoms) == 0:
        print(f"❌ Error: Could not identify ligand in the reference structure.")
        return np.array([])

    # Compute distances between protein and ligand atoms
    distances = md.compute_distances(ref_traj, np.array(np.meshgrid(protein_atoms, ligand_atoms)).T.reshape(-1, 2))
    distances = distances.reshape(len(protein_atoms), len(ligand_atoms))

    # Identify pocket residues within the cutoff distance
    pocket_atoms = protein_atoms[np.any(distances < cutoff, axis=1)]

    if len(pocket_atoms) == 0:
        print(f"⚠ Warning: No pocket residues found within {cutoff} nm.")

    return pocket_atoms

def main():
    parser = argparse.ArgumentParser(description="Compute RMSD between two PDB structures.")
    parser.add_argument("input_pdb", type=str, help="Path to the reference PDB file")
    parser.add_argument("target_pdb", type=str, help="Path to the target (aligned) PDB file")
    parser.add_argument("--cutoff", type=float, default=1.0, help="Cutoff distance for pocket residue identification (default: 1.0 nm)")
    args = parser.parse_args()

    print(f"📂 Loading reference PDB: {args.input_pdb}")
    try:
        ref_traj = md.load(args.input_pdb)
        if ref_traj.n_frames == 0:
            raise ValueError("Reference PDB is empty.")
    except Exception as e:
        print(f"❌ Error: Failed to load reference PDB ({args.input_pdb}): {e}")
        return

    print(f"📂 Loading target PDB: {args.target_pdb}")
    try:
        target_traj = md.load(args.target_pdb)
        if target_traj.n_frames == 0:
            raise ValueError("Target PDB is empty.")
    except Exception as e:
        print(f"❌ Error: Failed to load target PDB ({args.target_pdb}): {e}")
        return

    # Reference PDB: Get indices for different regions
    ref_backbone_idx = ref_traj.topology.select("protein and backbone and not element H")
    ref_ligand_idx = ref_traj.topology.select("not protein and not element H")
    #ref_ligand_idx = ref_traj.topology.select("not protein")
    pocket_atoms = identify_pocket_residues(ref_traj, args.cutoff)

    # Target PDB: Get indices for different regions
    target_backbone_idx = target_traj.topology.select("protein and backbone and not element H")
    target_ligand_idx = target_traj.topology.select("not protein and not element H")
    #target_ligand_idx = target_traj.topology.select("not protein")

    # Ensure selections are not empty before slicing
    description = "Protein Backbone"
    if len(ref_backbone_idx) > 0 and len(target_backbone_idx) > 0:
        ref_backbone = ref_traj.atom_slice(ref_backbone_idx)
        target_backbone = target_traj.atom_slice(target_backbone_idx)
        compute_rmsd(ref_backbone, target_backbone, description)
    else:
        print("⚠ Warning: Could not identify protein backbone in the target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")

    description = "Ligand"
    if len(ref_ligand_idx) > 0 and len(target_ligand_idx) > 0:
        ref_ligand = ref_traj.atom_slice(ref_ligand_idx)
        target_ligand = target_traj.atom_slice(target_ligand_idx)
        compute_rmsd(ref_ligand, target_ligand, description)
    else:
        print("⚠ Warning: Could not identify ligand in the target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")

    description = "Protein Pocket"
    if len(ref_backbone_idx) > 0 and len(pocket_atoms) > 0 and len(target_backbone_idx) > 0:
        ref_pocket = ref_traj.atom_slice(pocket_atoms)
        target_pocket = target_traj.atom_slice(pocket_atoms)
        compute_rmsd(ref_pocket, target_pocket, description)
    else:
        print("⚠ Warning: No protein pocket residues detected in the target PDB, skipping RMSD.")
        print(f"📊 {description} RMSD (Direct): None")
        print(f"📊 {description} RMSD (MDTraj): None")

if __name__ == "__main__":
    main()
