#!/usr/bin/env python3

import argparse
import sys

def extract_atom_info(line):
    return {
        "record": line[:6].strip(),
        "atom": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "coords": line[30:54]
    }

def atom_match(atom1, atom2):
    return atom1[0] == atom2[0]  # compare only the first character (element type)

def update_atom_blocks(pdb1_lines, pdb2_lines):
    atom_lines1 = [l for l in pdb1_lines if l.startswith(('ATOM', 'HETATM'))]
    atom_lines2 = [l for l in pdb2_lines if l.startswith(('ATOM', 'HETATM'))]

    if len(atom_lines1) != len(atom_lines2):
        sys.exit(f"❌ Atom count mismatch: {len(atom_lines1)} in pdb1 vs {len(atom_lines2)} in --ref template.")

    updated_atoms = []
    for i, (line1, line2) in enumerate(zip(atom_lines1, atom_lines2), start=1):
        info1 = extract_atom_info(line1)
        info2 = extract_atom_info(line2)

        if info1["resname"] != info2["resname"] or not atom_match(info1["atom"], info2["atom"]):
            print(f"❌ Mismatch at line {i}:")
            print(f"  pdb1: atom='{info1['atom']}' residue='{info1['resname']}'")
            print(f"  ref : atom='{info2['atom']}' residue='{info2['resname']}'")
            sys.exit("\n❌ Cannot update. Residue or atom type mismatch found.")

        # Replace coordinates from pdb1 into metadata from pdb2
        new_line = line2[:30] + info1["coords"] + line2[54:]
        updated_atoms.append(new_line)

    return updated_atoms

def merge_with_non_atoms(pdb2_lines, updated_atoms):
    output_lines = []
    atom_index = 0
    for line in pdb2_lines:
        if line.startswith(('ATOM', 'HETATM')):
            output_lines.append(updated_atoms[atom_index])
            atom_index += 1
        else:
            output_lines.append(line)
    return output_lines

def main():
    parser = argparse.ArgumentParser(
        description="Update coordinates in --ref PDB using pdb1. Match residues exactly, and atom types by first letter."
    )
    parser.add_argument("pdb1", help="PDB file with new coordinates (e.g., optimized structure)")
    parser.add_argument("--ref", required=True, help="Reference PDB file with atom metadata to preserve")
    parser.add_argument("-o", "--output", required=True, help="Output PDB file")
    args = parser.parse_args()

    with open(args.pdb1) as f1, open(args.ref) as f2:
        pdb1_lines = f1.readlines()
        pdb2_lines = f2.readlines()

    updated_atoms = update_atom_blocks(pdb1_lines, pdb2_lines)
    updated_lines = merge_with_non_atoms(pdb2_lines, updated_atoms)

    with open(args.output, 'w') as fout:
        for line in updated_lines:
            fout.write(line.rstrip() + "\n")

    print(f"✅ Coordinates from '{args.pdb1}' successfully applied to '{args.ref}' and saved to '{args.output}'.")

if __name__ == "__main__":
    main()

