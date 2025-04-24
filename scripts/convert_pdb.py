#!/usr/bin/env python3

import argparse
import sys

def get_coord_key(line):
    """Return coordinate block (X, Y, Z) as a string key."""
    return line[30:54]

def extract_atom_lines(path):
    """Extract ATOM/HETATM lines from a PDB file."""
    with open(path) as f:
        return [line for line in f if line.startswith(("ATOM", "HETATM"))]

def extract_serial(line):
    """Extract and return atom serial number (as int)."""
    return int(line[6:11])

def convert_pdb(pdb1_path, pdb2_path, pdb3_path, output_stream):
    pdb1_atoms = extract_atom_lines(pdb1_path)
    pdb2_atoms = extract_atom_lines(pdb2_path)
    pdb3_atoms = extract_atom_lines(pdb3_path)

    if not (len(pdb1_atoms) == len(pdb2_atoms) == len(pdb3_atoms)):
        raise ValueError("Fatal error: All three PDB files must have the same number of ATOM/HETATM lines.")

    # Step 1: Build coordinate → index map from PDB1
    coord_to_index = {get_coord_key(line): i for i, line in enumerate(pdb1_atoms)}

    # Step 2: Reorder PDB2 lines to match the atom order of PDB1
    reordered_pdb2 = [None] * len(pdb2_atoms)
    for i, line in enumerate(pdb2_atoms):
        coord = get_coord_key(line)
        if coord not in coord_to_index:
            raise ValueError(f"Fatal error: Coordinate {coord} in PDB2 not found in PDB1.")
        target_index = coord_to_index[coord]
        reordered_pdb2[target_index] = line

    # Step 3: Build the new PDB4 lines with atom identity check
    pdb4_lines = []

    for i, line3 in enumerate(pdb3_atoms):
        line2 = reordered_pdb2[i]

        # Extract atom name and residue name
        atom_name_2 = line2[12:16].strip()
        res_name_2 = line2[17:20].strip()
        atom_name_3 = line3[12:16].strip()
        res_name_3 = line3[17:20].strip()

        if atom_name_2 != atom_name_3 or res_name_2 != res_name_3:
            raise ValueError(f"❌ Mismatch at atom {i}:\n"
                             f"    PDB2: {atom_name_2} {res_name_2}\n"
                             f"    PDB3: {atom_name_3} {res_name_3}")

        # Build new line: keep atom identity from line2, use coordinates from line3
        new_line = line2[:30] + line3[30:54] + line2[54:]
        pdb4_lines.append(new_line)

    # Step 4: Sort lines by atom serial number
    sorted_pdb4 = sorted(pdb4_lines, key=extract_serial)

    # Output
    for line in sorted_pdb4:
        output_stream.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Generate pdb4 by combining coordinates from pdb3 and atom info from pdb2 by following pdb1-pdb2 coordinate mapping," +
                    "like we know pdb2 <- pdb1, then we need to know pdb4 <- pdb3 with the same coordinate mapping",
        usage="python convert_pdb.py pdb1_native_modified.pdb pdb2_native_ori.pdb pdb3.pdb > pdb4.pdb"
    )
    parser.add_argument("pdb1", help="PDB with different atom orders with PDB2, e.g., native_reordered.pdb")
    parser.add_argument("pdb2", help="PDB with desired atom info (existing a coords mapping with PDB1, e.g., native.pdb)")
    parser.add_argument("pdb3", help="PDB with desired coordinates and the same atom type & residue name with PDB2")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    convert_pdb(args.pdb1, args.pdb2, args.pdb3, sys.stdout)

if __name__ == "__main__":
    main()

