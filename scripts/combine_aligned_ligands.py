import os
import argparse
import mdtraj as md

def combine_aligned_ligands(root_dir, output_filename):
    """
    Combines all `aligned_ligand.pdb` files into a single PDB file.

    Parameters:
    - root_dir (str): The root directory containing `aligned_ligand.pdb` files.
    - output_filename (str): The final output file name for the combined ligands.
    """
    ligand_trajs = []

    # Walk through all directories and find aligned_ligand.pdb files
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename == "aligned_ligand.pdb":  # Find all ligand PDBs
                ligand_pdb_path = os.path.join(dirpath, filename)
                try:
                    traj = md.load(ligand_pdb_path)
                    ligand_trajs.append(traj)
                    print(f"📂 Loaded: {ligand_pdb_path}")
                except Exception as e:
                    print(f"❌ Error loading {ligand_pdb_path}: {e}")

    if not ligand_trajs:
        print("⚠ No `aligned_ligand.pdb` files found. Exiting.")
        return

    # Concatenate all ligand trajectories into one
    combined_traj = md.join(ligand_trajs)

    # Save the combined ligands as a new PDB file
    combined_traj.save(output_filename)
    print(f"✅ Combined ligand PDB saved as: {output_filename}")

def main():
    """Main function to handle argument parsing."""
    parser = argparse.ArgumentParser(
        description="Combines all `aligned_ligand.pdb` files from a directory into a single PDB file."
    )
    parser.add_argument(
        "--dir", default=".", help="Root directory that has the folders containing 'aligned_ligand.pdb' files (default: current directory)."
    )
    parser.add_argument(
        "-o", "--output", default="combined_aligned_ligands.pdb", 
        help="Output PDB file name (default: 'combined_aligned_ligands.pdb')."
    )

    args = parser.parse_args()

    # Run the function with the provided arguments
    combine_aligned_ligands(args.dir, args.output)

if __name__ == "__main__":
    main()

