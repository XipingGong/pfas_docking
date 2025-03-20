import sys
import argparse
import mdtraj as md

# Define expected heavy atom count for each amino acid (excluding hydrogen)
HEAVY_ATOM_COUNT = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9,
    "GLY": 4, "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11,
    "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7
}

def check_pdb_ligands_and_chains(pdb_file, ligand_id=None):
    """
    Analyze a PDB file to identify protein and ligand chains, including water molecules, printing relevant details.

    :param pdb_file: Path to the PDB file.
    :param ligand_id: Specific ligand ID to check for multiple presence (optional).
    :return: Prints details of protein chains and all detected ligands, including water and missing heavy atoms.
    """
    # Load PDB file
    traj = md.load(pdb_file)
    print(f"{pdb_file}: {traj}")

    # Identify protein chains and count residues and atoms
    protein_info = {}
    missing_atoms = {}
    for residue in traj.top.residues:
        if residue.is_protein:
            chain_id = residue.chain.chain_id
            resname = residue.name.upper()
            resid = residue.index  # Residue index
            resSeq = residue.resSeq  # Residue sequence number
            
            if chain_id not in protein_info:
                protein_info[chain_id] = {"residues": 0, "atoms": 0}
                missing_atoms[chain_id] = []  # Track missing atoms
            protein_info[chain_id]["residues"] += 1
            protein_info[chain_id]["atoms"] += len(list(residue.atoms))
            
            # Check for missing heavy atoms
            expected_atom_count = HEAVY_ATOM_COUNT.get(resname, 0)
            actual_atom_count = len(list(residue.atoms))
            if actual_atom_count < expected_atom_count:
                missing_atoms[chain_id].append((resid, resSeq, expected_atom_count - actual_atom_count))

    if protein_info:
        print("\nProtein Info:")
        for chain_id, info in sorted(protein_info.items()):
            print(f"  - Chain '{chain_id}':")
            print(f"    - Number of residues: {info['residues']}")
            print(f"    - Number of atoms: {info['atoms']}")
            if missing_atoms[chain_id]:
                print(f"    - {pdb_file}: Missing heavy atoms: {len(missing_atoms[chain_id]) > 0}")
                for resid, resSeq, missing_count in missing_atoms[chain_id]:
                    print(f"      - Resid: {resid}, ResSeq: {resSeq}, Missing Atoms: {missing_count}")
        print(f"  - {pdb_file}: Multiple protein chains present: {len(protein_info) > 1}")
    else:
        print("No protein chains found in the PDB file.")

    # Identify unique ligand residues (non-protein, including water)
    ligands = {}
    for residue in traj.top.residues:
        if not residue.is_protein:
            chain_id = residue.chain.chain_id
            if residue.name not in ligands:
                ligands[residue.name] = {"chains": set(), "residues": 0, "atoms": 0}
            ligands[residue.name]["chains"].add(chain_id)
            ligands[residue.name]["residues"] += 1
            ligands[residue.name]["atoms"] += len(list(residue.atoms))

    if ligands:
        print("\nLigand Info (including water):")
        for ligand_name, info in sorted(ligands.items()):
            print(f"  - Ligand '{ligand_name}':")
            print(f"    - Chain IDs: {sorted(info['chains'])}")
            print(f"    - Number of residues: {info['residues']}")
            print(f"    - Number of atoms: {info['atoms']}")
            if ligand_name == ligand_id:
                print(f"    - {pdb_file}: Multiple '{ligand_id}' ligand residues present: {info['residues'] > 1}")
    else:
        print("No ligands found in the PDB file.")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze PDB file for protein and ligand information.")
    parser.add_argument("pdb_file", type=str, help="Path to the PDB file.")
    parser.add_argument("--ligand_id", type=str, default='8PF', help="Ligand ID to check for multiple presence.")
    
    args = parser.parse_args()
    check_pdb_ligands_and_chains(args.pdb_file, args.ligand_id)

