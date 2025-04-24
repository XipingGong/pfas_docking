import re
import argparse

def convert_ligand_file(input_file):
    """
    Converts a ligand file from "Ligand ID: X ; PDBID: Y ;" format 
    to "Y --ligand_id X" and prints the results.

    Parameters:
    ----------
    input_file : str
        Path to the input text file.
    """
    with open(input_file, "r") as infile:
        for line in infile:
            ligand_match = re.search(r'Ligand ID:\s*(\S+)\s*;', line)
            pdb_match = re.search(r'PDBID:\s*([\S ]+)\s*;', line)

            if ligand_match and pdb_match:
                ligand_id = ligand_match.group(1)
                pdb_ids = pdb_match.group(1).split()

                # Ignore if PDB ID is "not_available"
                for pdb_id in pdb_ids:
                    if pdb_id.lower() != "not_available":
                        print(f"{pdb_id}_{ligand_id}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ligand file format and print the results.")
    parser.add_argument("input_file", type=str, help="Path to the input text file.")

    args = parser.parse_args()
    convert_ligand_file(args.input_file)

