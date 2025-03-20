#!/usr/bin/env python3

import re
import sys
import argparse
from rdkit import Chem

def parse_and_filter_ligands(file_path, min_fluorine):
    """
    Extracts ligand ID and SMILES from a CIF file and filters ligands based on a minimum number of fluorine (F) atoms.

    Parameters:
        file_path (str): Path to the CIF file.
        min_fluorine (int): Minimum number of fluorine atoms required in a ligand.

    Returns:
        list: A sorted list of tuples containing (Ligand ID, SMILES, Fluorine Count) for ligands meeting the fluorine threshold.
    """
    id_pattern = re.compile(r"_chem_comp\.id\s+(\S+)")
    smiles_pattern = re.compile(r"^\S+\s+SMILES_CANONICAL\s+CACTVS\s+\S+\s+\"([^\"]+)\"")  # Prioritize CACTVS

    filtered_ligands = []
    ligand_id, smiles = None, None

    with open(file_path, "r") as f:
        for line in f:
            id_match = id_pattern.match(line)
            smiles_match = smiles_pattern.search(line)

            if id_match:
                ligand_id = id_match.group(1)
            if smiles_match:
                smiles = smiles_match.group(1)

            if ligand_id and smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        f_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "F")
                        if f_count >= min_fluorine:
                            filtered_ligands.append((ligand_id, smiles, f_count))
                except Exception as e:
                    print(f"⚠️ Skipping invalid SMILES for Ligand ID {ligand_id}: {smiles} - Error: {e}")
                ligand_id, smiles = None, None  # Reset for next ligand

    # Sort ligands by fluorine count in descending order
    filtered_ligands.sort(key=lambda x: x[2], reverse=True)
    return filtered_ligands

def main():
    """
    Main execution function.

    Usage:
        python get_fluorine_ligand.py <components.cif> <min_fluorine_count>

    Arguments:
        components.cif : The input CIF file containing ligand structures. This file can be downloaded from the Protein Data Bank (PDB), typically containing metadata for small molecules, ligands, and their corresponding chemical structures. It follows the Crystallographic Information File (CIF) format, which is used for describing molecular and crystal structures.

        The components.cif file can be downloaded from the PDB Chemical Component Dictionary:
        https://files.wwpdb.org/pub/pdb/data/monomers/components.cif

        min_fluorine_count : The minimum number of fluorine atoms required in a ligand to be included in the output. Only ligands with at least this number of fluorine atoms will be considered.

    Output:
        Prints the Ligand ID, corresponding SMILES, and the count of fluorine atoms for ligands meeting the fluorine threshold.
        The output is sorted in descending order by fluorine count.
    """
    parser = argparse.ArgumentParser(description="Extract ligand IDs and SMILES from a CIF file, filtering by fluorine count and sorting by highest fluorine content.")
    parser.add_argument("cif_file", type=str, help="Path to the CIF file containing ligand structures from the Protein Data Bank (PDB). This file can be downloaded from: https://files.wwpdb.org/pub/pdb/data/monomers/components.cif")
    parser.add_argument("min_fluorine_count", type=int, help="Minimum number of fluorine atoms required in a ligand.")
    args = parser.parse_args()

    filtered_ligands = parse_and_filter_ligands(args.cif_file, args.min_fluorine_count)

    if filtered_ligands:
        print(f"# Ligand IDs - Found {len(filtered_ligands)} ligands with at least {args.min_fluorine_count} fluorine atoms:")
        for lig_id, smiles, f_count in filtered_ligands:
            print(f"Ligand ID: {lig_id} ; SMILES: {smiles} ; Fluorine Count: {f_count} ;")
    else:
        print(f"# No ligands found with at least {args.min_fluorine_count} fluorine atoms.")

if __name__ == "__main__":
    main()

