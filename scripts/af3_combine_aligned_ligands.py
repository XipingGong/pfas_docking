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

    # Collect all aligned_ligand.pdb paths
    aligned_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "aligned_ligand.pdb" in filenames:
            aligned_files.append(os.path.join(dirpath, "aligned_ligand.pdb"))

    # Sort with "best_pose" first
    aligned_files.sort(key=lambda x: (0 if "best_pose" in x else 1, x))

    # Load all sorted aligned_ligand.pdb files
    ligand_trajs = []
    for path in aligned_files:
        try:
            ligand_trajs.append(md.load(path))
            print(f"📂 Loaded: {path}")
        except Exception as e:
            print(f"❌ Error loading {path}: {e}")

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

