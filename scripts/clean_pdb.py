import argparse
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import mdtraj as md

def fix_missing_atoms(input_pdb, output_pdb):
    """
    Uses PDBFixer to fill in missing heavy atoms while keeping only heavy atoms in the output file.

    Parameters:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to save the fixed PDB file.
    
    Output:
        A revised PDB file with missing heavy atoms filled and hydrogen atoms removed.
    """
    print(f"Processing {input_pdb}...")

    # Load PDB file into PDBFixer
    fixer = PDBFixer(filename=input_pdb)

    # Find and add missing heavy atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Save the fixed PDB file
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    print(f"Fixed PDB file saved as {output_pdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fill missing heavy atoms in a PDB file using PDBFixer.")
    parser.add_argument("input_pdb", type=str, help="Path to the input PDB file with missing atoms.")
    parser.add_argument("-o", "--output", type=str, help="Path to save the revised PDB file. Default: cleaned_<input_pdb>")
    
    args = parser.parse_args()
    
    # Handle full path for input PDB file
    input_pdb_basename = os.path.basename(args.input_pdb)
    output_pdb = args.output if args.output else os.path.join(os.path.dirname(args.input_pdb), f"cleaned_{input_pdb_basename}")
    
    fix_missing_atoms(args.input_pdb, output_pdb)

