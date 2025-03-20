#!/usr/bin/env python3 

import re
import argparse
import pandas as pd
from rdkit import Chem

def filter_ligands(log_file, smarts_str):
    """
    Reads a .log file to parse ligands and SMILES, filters for those
    containing fully fluorinated methyl/methylene units based on the
    specified SMARTS pattern.

    Parameters
    ----------
    log_file : str
        Path to the log file that contains lines with ligand ID and SMILES.
    smarts_str : str
        A SMARTS pattern used to detect the desired fully fluorinated units.
    """
    # Read the file lines
    with open(log_file, "r") as file:
        data = file.readlines()

    # Extract Ligand ID, SMILES, and PDB IDs using regex
    ligand_dict = {}
    ligand_regex = re.compile(r'Ligand ID:\s*(\S+)\s*;\s*SMILES:\s*([^;]+)')
    pdb_regex = re.compile(r'PDBID:\s*([^;]+)')

    current_ligand = None

    for line in data:
        ligand_match = ligand_regex.search(line)
        if ligand_match:
            ligand_id = ligand_match.group(1).strip()
            smiles = ligand_match.group(2).strip()
            if ligand_id not in ligand_dict:
                ligand_dict[ligand_id] = {"SMILES": smiles, "PDB IDs": set()}
            current_ligand = ligand_id

        pdb_match = pdb_regex.search(line)
        if pdb_match and current_ligand:
            pdb_ids = pdb_match.group(1).strip().split()
            ligand_dict[current_ligand]["PDB IDs"].update(pdb_ids)

    if not ligand_dict:
        print("No valid ligand data found in the log file.")
        return

    # Convert to DataFrame
    ligands = [(lig_id, info["SMILES"], " ".join(sorted(info["PDB IDs"]))) 
               for lig_id, info in ligand_dict.items()]
    df = pd.DataFrame(ligands, columns=["Ligand ID", "SMILES", "PDB IDs"])

    # Convert SMILES to RDKit Molecule objects
    def safe_mol_from_smiles(smile):
        """Helper function to convert SMILES safely."""
        mol = Chem.MolFromSmiles(smile)
        return mol if mol else None

    df["Molecule"] = df["SMILES"].apply(safe_mol_from_smiles)

    # Drop invalid molecules
    df = df.dropna(subset=["Molecule"])

    # Convert the SMARTS string into an RDKit pattern
    pfas_pattern = Chem.MolFromSmarts(smarts_str)
    if pfas_pattern is None:
        raise ValueError("Invalid SMARTS pattern. Check syntax.")

    # Check if the molecule has the specified substructure
    df["Has_Fully_Fluorinated_Unit"] = df["Molecule"].apply(lambda mol: mol.HasSubstructMatch(pfas_pattern))
    df_fully_fluorinated = df[df["Has_Fully_Fluorinated_Unit"] == True]

    # Print summary
    print(f"# Number of ligands matching the pattern {smarts_str}: {len(df_fully_fluorinated)} (Total = {len(df)})")
    for _, row in df_fully_fluorinated.iterrows():
        print(f"Ligand ID: {row['Ligand ID']} ; PDBID: {row['PDB IDs']} ;")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter ligands from a log file based on a SMARTS pattern.",
                                     epilog="Example usage: python filter_ligands.py path/to/log_file.log '[CX4]([!#1])(F)(F)-[CX4](F)([!#1])([!#1])' # for R-CF2-CF(R')(R\") where R, R', R\" ≠ H")
    parser.add_argument("log_file", type=str, help="Path to the log file containing ligand data.")
    parser.add_argument("smarts_str", type=str, help="SMARTS pattern to filter ligands.")

    args = parser.parse_args()
    filter_ligands(args.log_file, args.smarts_str)

