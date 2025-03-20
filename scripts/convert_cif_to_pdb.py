import os
import argparse
import subprocess

def convert_cif_to_pdb(root_dir):
    """
    Finds all 'model.cif' files in subdirectories and converts them to 'model.pdb' using Open Babel (obabel).

    Parameters:
    - root_dir (str): The root directory where the search begins.
    """
    if not os.path.exists(root_dir):
        print(f"❌ Error: The specified directory does not exist: {root_dir}")
        return

    converted_count = 0

    # Walk through all directories and subdirectories
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower() == "model.cif":  # Find all model.cif files
                cif_path = os.path.join(dirpath, filename)  # Full path to model.cif
                pdb_path = os.path.join(dirpath, "model.pdb")  # Output model.pdb path

                # Run Open Babel (obabel) to convert CIF to PDB
                command = ["obabel", cif_path, "-O", pdb_path]
                try:
                    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    print(f"✅ Converted: {cif_path} → {pdb_path}")
                    converted_count += 1
                except subprocess.CalledProcessError as e:
                    print(f"❌ Error: Failed to convert {cif_path}\nError: {e}")

    if converted_count == 0:
        print("⚠ No 'model.cif' files found for conversion.")

def main():
    """Main function to handle argument parsing."""
    parser = argparse.ArgumentParser(
        description="Convert all 'model.cif' files in a directory to 'model.pdb' using Open Babel."
    )
    parser.add_argument(
        "root_dir", nargs="?", default=".",
        help="Root directory to search for 'model.cif' files (default: current directory)."
    )

    args = parser.parse_args()

    # Run the conversion function
    convert_cif_to_pdb(args.root_dir)

if __name__ == "__main__":
    main()

