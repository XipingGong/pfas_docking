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

def get_atom_type(line):
    """Extract atom type (element symbol) from column 77-78."""
    return line[76:78].strip()

def atom_types_match(line1, line2):
    """Check if the first letter of atom types match."""
    return get_atom_type(line1)[:1] == get_atom_type(line2)[:1]

def merge_atom_info(coord_line, info_line):
    """Replace the atom info in info_line with coordinates from coord_line."""
    return info_line[:30] + coord_line[30:54] + info_line[54:]

def reorder_pdb_by_coord_mapping(pdb1_path, pdb2_path, pdb3_path, output_stream, pdbinfo_path="USE_PDB2"):
    pdb1_atoms = extract_atom_lines(pdb1_path)
    pdb2_atoms = extract_atom_lines(pdb2_path)
    pdb3_atoms = extract_atom_lines(pdb3_path)

    if not (len(pdb1_atoms) == len(pdb2_atoms) == len(pdb3_atoms)):
        raise ValueError("❌ All input PDBs must have the same number of ATOM/HETATM lines.")

    # Step 1: Build mapping from coordinates in PDB1 → index
    coord_to_index_pdb1 = {get_coord_key(line): j for j, line in enumerate(pdb1_atoms)}

    # Step 2: Build reorder map from PDB1 → PDB2
    reorder_map = []
    for i, line in enumerate(pdb2_atoms):
        coord = get_coord_key(line)
        if coord not in coord_to_index_pdb1:
            raise ValueError(f"❌ Coordinate {coord.strip()} from PDB2 not found in PDB1.")
        reorder_map.append(coord_to_index_pdb1[coord])

    # Step 3: Determine which atom info source to use
    info_atoms = pdb2_atoms if pdbinfo_path == "USE_PDB2" else extract_atom_lines(pdbinfo_path)

    if len(info_atoms) != len(pdb3_atoms):
        raise ValueError("❌ PDBINFO must have the same number of atoms as PDB3.")

    # Step 4: Merge atom info with reordered coordinates
    final_lines = []
    for i, j in enumerate(reorder_map):
        coord_line = pdb3_atoms[j]
        info_line = info_atoms[i]
        if atom_types_match(coord_line, info_line):
            final_lines.append(merge_atom_info(coord_line, info_line))
        else:
            raise ValueError(
                f"❌ Atom type mismatch at index {i}: {get_atom_type(coord_line)} vs {get_atom_type(info_line)}"
            )

    # Step 5: Write output
    for line in final_lines:
        output_stream.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Reorder PDB3 → PDB4 using mapping from PDB1 → PDB2 (based on coordinate match). Optionally replace atom info.",
        usage="python reorder_pdb_by_coord_mapping pdb1.pdb pdb2.pdb pdb3.pdb [--pdbinfo ref.pdb] > pdb4.pdb"
    )
    parser.add_argument("pdb1", help="Reference PDB defining desired atom order")
    parser.add_argument("pdb2", help="PDB with same coordinates but different atom order")
    parser.add_argument("pdb3", help="PDB whose coordinates will be reordered to match PDB1")
    parser.add_argument(
        "--pdbinfo",
        default="USE_PDB2",
        help="Optional PDB file to take atom info (excluding coordinates). Default is to use pdb2."
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    reorder_pdb_by_coord_mapping(
        pdb1_path=args.pdb1,
        pdb2_path=args.pdb2,
        pdb3_path=args.pdb3,
        output_stream=sys.stdout,
        pdbinfo_path=args.pdbinfo
    )

if __name__ == "__main__":
    main()

