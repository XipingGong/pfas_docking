#!/usr/bin/env python3
import argparse
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import SanitizeMol, SanitizeFlags

def read_pdb_blocks(path):
    header, atoms, footer = [], [], []
    mode = "header"
    with open(path) as f:
        for line in f:
            if mode == "header":
                if line.startswith(("ATOM  ", "HETATM")):
                    mode = "atoms"
                    atoms.append(line)
                else:
                    header.append(line)
            elif mode == "atoms":
                if line.startswith(("ATOM  ", "HETATM")):
                    atoms.append(line)
                else:
                    mode = "footer"
                    footer.append(line)
            else:
                footer.append(line)
    return header, atoms, footer

def parse_xyz(line):
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    return x, y, z

def main():
    parser = argparse.ArgumentParser(
        description="Rename REF atom-info + splice MOV coords → stdout PDB"
    )
    parser.add_argument(
        "--ref", required=True,
        help="Reference PDB (native_ligand.pdb) providing atom info"
    )
    parser.add_argument(
        "mov",
        help="Moving PDB (lig.pdb) providing coordinates"
    )
    args = parser.parse_args()

    # Load molecules in RDKit without valence checking
    ref_mol = Chem.MolFromPDBFile(args.ref, removeHs=False, sanitize=False)
    mov_mol = Chem.MolFromPDBFile(args.mov, removeHs=False, sanitize=False)
    if ref_mol is None or mov_mol is None:
        sys.exit("ERROR: RDKit failed to parse one of the PDBs")

    # Sanitize without the valence/properties check
    flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES
    SanitizeMol(ref_mol, flags)
    SanitizeMol(mov_mol, flags)

    # Get atom-index mapping: ref_idx -> mov_idx
    mapping = mov_mol.GetSubstructMatch(ref_mol)
    if len(mapping) != ref_mol.GetNumAtoms():
        sys.exit("ERROR: mapping length mismatch between ref and mov PDBs")

    # Read PDB text blocks
    ref_header, ref_atoms, ref_footer = read_pdb_blocks(args.ref)
    _, mov_atoms, _ = read_pdb_blocks(args.mov)
    if len(ref_atoms) != len(mov_atoms):
        sys.exit("ERROR: atom-line counts differ between ref and mov PDBs")

    # Extract moving coords
    mov_xyz = np.array([parse_xyz(l) for l in mov_atoms])

    # Splice coords into reference lines
    for i, line in enumerate(ref_atoms):
        j = mapping[i]
        x, y, z = mov_xyz[j]
        new_line = (
            line[:30]
            + f"{x:8.3f}{y:8.3f}{z:8.3f}"
            + line[54:]
        )
        sys.stdout.write(new_line)

    # Write footer
    sys.stdout.writelines(ref_footer)

if __name__ == "__main__":
    main()

