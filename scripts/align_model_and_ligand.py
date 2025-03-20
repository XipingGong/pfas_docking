import os
import argparse
import mdtraj as md
import numpy as np

def align_model_and_extract_ligand(pdb_path, ref_pdb, cutoff=1.0):
    """
    Aligns the protein model using pocket residues from a reference PDB and extracts the ligand.

    Parameters:
    - pdb_path (str): Path to the input model.pdb file.
    - ref_pdb (str): Path to the reference PDB file.
    - cutoff (float): Distance cutoff in nanometers (default: 1.0 nm = 10 Å).

    Outputs:
    - Saves `aligned_model.pdb`
    - Saves `aligned_ligand.pdb`
    """

    # Load the reference PDB and target PDB
    print(f"📂 Loading Reference PDB: {ref_pdb}")
    ref_traj = md.load(ref_pdb)
    print(ref_traj)  # Print summary of the reference trajectory

    print(f"📂 Loading Target PDB: {pdb_path}")
    traj = md.load(pdb_path)
    print(traj)  # Print summary of the target trajectory


    # Select protein backbone (chain A) and ligand (chain B)
    protein_atoms = ref_traj.topology.select("protein and backbone and not element H")
    ligand_atoms = ref_traj.topology.select("not protein and not element H")  # Ligand (chain B)

    if len(protein_atoms) == 0 or len(ligand_atoms) == 0:
        print(f"❌ Error: Could not identify both protein and ligand in {pdb_path}")
        return

    # Compute distances between protein atoms and ligand atoms
    pairs = np.array([[p, l] for p in protein_atoms for l in ligand_atoms])
    distances = md.compute_distances(ref_traj, pairs).reshape(len(protein_atoms), len(ligand_atoms))

    # Identify pocket residues (within cutoff distance)
    pocket_mask = np.any(distances < cutoff, axis=1)
    pocket_atoms = protein_atoms[pocket_mask]

    if len(pocket_atoms) == 0:
        print(f"⚠ Warning: No pocket residues found for {pdb_path}. Alignment may be inaccurate.")

    # Align the model to the reference using pocket residues
    traj.superpose(ref_traj, atom_indices=pocket_atoms)

    # Define output filenames
    output_dir = os.path.dirname(pdb_path)
    aligned_model_pdb = os.path.join(output_dir, "aligned_model.pdb")
    aligned_ligand_pdb = os.path.join(output_dir, "aligned_ligand.pdb")

    # Save the aligned full structure
    traj.save(aligned_model_pdb)
    print(f"✅ Saved aligned model: {aligned_model_pdb}")

    # Compute RMSD of protein
    idx = traj.topology.select(f"protein and backbone and not element H")
    target_protein_traj = traj.atom_slice(idx)
    idx = ref_traj.topology.select(f"protein and backbone and not element H")
    ref_protein_traj = ref_traj.atom_slice(idx)
    rmsd_values = md.rmsd(target_protein_traj, ref_protein_traj)
    print(f"📊 Pocket-Aligned Protein Backbone RMSD (MDTraj): {rmsd_values} nm")

    # Compute RMSD of protein pocket
    target_pocket_traj = traj.atom_slice(pocket_atoms)
    ref_pocket_traj = ref_traj.atom_slice(pocket_atoms)
    rmsd_values = md.rmsd(target_pocket_traj, ref_pocket_traj)
    print(f"📊 Pocket-Aligned Protein Pocket RMSD (MDTraj): {rmsd_values} nm")

    # Extract only the ligand and save separately
    ligand_traj = traj.atom_slice(ligand_atoms)
    ligand_traj.save(aligned_ligand_pdb)
    print(f"✅ Saved aligned ligand: {aligned_ligand_pdb}")

def process_all_models(root_dir, ref_pdb):
    """
    Searches for all `model.pdb` files and processes them using the given reference PDB.
    """
    if not os.path.exists(ref_pdb):
        print(f"❌ Error: Reference PDB file not found: {ref_pdb}")
        return

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename == "model.pdb":  # Process only model.pdb files
                pdb_path = os.path.join(dirpath, filename)
                align_model_and_extract_ligand(pdb_path, ref_pdb)

def main():
    """Main function to handle argument parsing."""
    parser = argparse.ArgumentParser(
        description="Aligns model PDBs to a reference PDB using pocket residues and extracts ligands."
    )
    parser.add_argument(
        "--dir", nargs="?", default=".",
        help="Root directory containing folders that have the 'model.pdb' files (default: current directory)."
    )
    parser.add_argument(
        "ref_pdb", help="Path to the reference PDB file ('input.pdb')."
    )

    args = parser.parse_args()

    # Run the process with the correct argument
    process_all_models(args.dir, args.ref_pdb)

if __name__ == "__main__":
    main()

