#!/usr/bin/env python3

import argparse

def merge_pdb_lines(pdb1, pdb2, output):
    with open(pdb1, 'r') as f1, open(pdb2, 'r') as f2, open(output, 'w') as fout:
        for line in f1:
            if line.startswith(("ATOM", "HETATM")):
                fout.write(line)
        for line in f2:
            if line.startswith(("ATOM", "HETATM")):
                fout.write(line)
        fout.write("END\n")
    print(f"✅ Merged PDB (with exact atom ordering) saved to {output}")

def main():
    parser = argparse.ArgumentParser(description="Safely merge two PDB files without changing atom order or names.")
    parser.add_argument("pdb1", help="First PDB file")
    parser.add_argument("pdb2", help="Second PDB file")
    parser.add_argument("-o", "--output", required=True, help="Output merged PDB file")
    args = parser.parse_args()

    merge_pdb_lines(args.pdb1, args.pdb2, args.output)

if __name__ == "__main__":
    main()

