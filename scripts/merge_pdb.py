#!/usr/bin/env python3

import argparse
import sys

def parse_pdb_lines(filepath):
    with open(filepath) as f:
        return f.readlines()

def is_atom_line(line):
    return line.startswith("ATOM") or line.startswith("HETATM")

def is_hydrogen(line):
    return line[76:78].strip() == 'H'

def format_atom_name(name):
    return name.rjust(4)  # Align to 4-character PDB convention

def get_element_symbol(line):
    return line[76:78].strip()

def get_atom_name(line):
    return line[12:16].strip()

def merge_pdb_lines(ref_lines, coord_lines):
    pdb3_lines = []
    ref_atoms = [line for line in ref_lines if is_atom_line(line) and not is_hydrogen(line)]
    ref_index = 0
    atom_serial = 1
    hydrogen_counter = 0

    hydrogen_names = [f"H{i}" if i > 0 else "H" for i in range(100)]

    for coord_line in coord_lines:
        if not is_atom_line(coord_line):
            pdb3_lines.append(coord_line.rstrip())
            continue

        is_H = is_hydrogen(coord_line)

        if not is_H:
            if ref_index >= len(ref_atoms):
                sys.stderr.write("⚠️ Warning: Not enough reference atoms to match non-H atoms in coord file.\n")
                break

            ref_line = ref_atoms[ref_index]
            ref_index += 1

            # Atom name consistency check (first letter of atom name)
            ref_elem = get_atom_name(ref_line)[0]
            coord_elem = get_atom_name(coord_line)[0]
            if ref_elem != coord_elem:
                sys.stderr.write(
                    f"❌ Error: Atom name mismatch at index {ref_index}: "
                    f"'{get_atom_name(ref_line)}' vs '{get_atom_name(coord_line)}'\n"
                )
                sys.exit(1)

            # Coordinate check (X, Y, Z must match exactly)
            if ref_line[30:54] != coord_line[30:54]:
                sys.stderr.write(
                    f"❌ Error: Coordinate mismatch at atom {ref_index}:\n"
                    f"    Ref:   {ref_line[30:54]}\n"
                    f"    Coord: {coord_line[30:54]}\n"
                )
                sys.exit(1)

            # Merge line with updated serial
            new_line = (
                ref_line[:6] +
                f"{atom_serial:5d} " +
                ref_line[12:16] +
                ref_line[16:30] +
                coord_line[30:54] +
                coord_line[54:76] +
                coord_line[76:]
            )
        else:
            # Generate hydrogen atom name (H, H1, H2, ...)
            h_name = hydrogen_names[hydrogen_counter] if hydrogen_counter < len(hydrogen_names) else f"H{hydrogen_counter}"
            atom_name = format_atom_name(h_name)
            hydrogen_counter += 1

            ref_base = ref_atoms[ref_index - 1] if ref_index > 0 else coord_line

            new_line = (
                ref_base[:6] +
                f"{atom_serial:5d} " +
                atom_name +
                ref_base[16:30] +
                coord_line[30:54] +
                coord_line[54:76] +
                coord_line[76:]
            )

        pdb3_lines.append(new_line.rstrip())
        atom_serial += 1

    return pdb3_lines

def main():
    parser = argparse.ArgumentParser(
        description="""Merge two PDB files:
- Uses pdb2 as the coordinate file (includes hydrogens).
- Uses metadata from pdb1 for non-H atoms.
- Ensures atom names start with the same element letter.
- Ensures matching coordinates for non-H atoms.
- Renames H atoms as H, H1, H2, ..., H999.
- Keeps all other lines (REMARK, TER, etc.) from pdb2.
- Renumbers atom serials sequentially.""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pdb1", help="Reference PDB (no H atoms, clean metadata)")
    parser.add_argument("pdb2", help="Coordinate PDB (includes H atoms)")

    args = parser.parse_args()

    ref_lines = parse_pdb_lines(args.pdb1)
    coord_lines = parse_pdb_lines(args.pdb2)

    merged = merge_pdb_lines(ref_lines, coord_lines)

    for line in merged:
        print(line)
    print("END")

if __name__ == "__main__":
    main()

