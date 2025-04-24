#!/usr/bin/env python3

import re
import argparse
from rdkit import Chem

def print_matching_lines(log_file, smarts_str):
    """
    Reads a log file and prints lines containing SMILES strings matching the SMARTS pattern.
    """
    ligand_regex = re.compile(r'SMILES:\s*([^;]+)')

    # Compile SMARTS
    pattern = Chem.MolFromSmarts(smarts_str)
    if pattern is None:
        raise ValueError("Invalid SMARTS pattern.")

    with open(log_file, 'r') as f:
        for line in f:
            match = ligand_regex.search(line)
            if match:
                smiles = match.group(1).strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(pattern):
                    print(line.strip())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Print lines from a log file where the SMILES match the given SMARTS pattern."
    )
    parser.add_argument("log_file", type=str, help="Path to the log file.")
    parser.add_argument("--pattern", required=True, type=str, help="SMARTS pattern to filter lines.")

    args = parser.parse_args()
    print_matching_lines(args.log_file, args.pattern)

