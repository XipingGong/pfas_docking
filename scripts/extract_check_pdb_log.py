import re
import argparse

def analyze_pdb_output(input_file):
    """
    Parses the log file and extracts information on missing atoms,
    multiple ligands, multiple protein chains, total number of atoms, and total number of residues for each PDB file.

    Parameters:
        input_file (str): Path to the log file.

    Prints:
        The PDB filename followed by:
        - Total number of atoms
        - Total number of residues
        - Missing heavy atoms (True/False)
        - Multiple ligand residues present (True/False)
        - Multiple protein chains present (True/False)
        - PDB type classification
    """
    with open(input_file, "r") as file:
        lines = file.readlines()
    
    results = []
    current_pdb = None
    ligand_id = None
    total_atoms = 0
    total_residues = 0
    missing_atoms = False
    multiple_ligands = False
    multiple_protein_chains = False

    for line in lines:
        # Identify the start of a new PDB section
        match = re.search(r'check_pdb\.py.*?(\S+\.pdb)\s+--ligand_id\s+(\S+)', line)
        if match:
            # Save the previous PDB results before starting a new one
            if current_pdb:
                pdb_type = classify_pdb_type(total_atoms, total_residues, multiple_ligands, multiple_protein_chains)
                results.append((current_pdb, total_atoms, total_residues, missing_atoms, multiple_ligands, multiple_protein_chains, pdb_type))
            
            # Reset values for the new PDB entry
            current_pdb = match.group(1)
            ligand_id = match.group(2)
            total_atoms = 0
            total_residues = 0
            missing_atoms = False
            multiple_ligands = False
            multiple_protein_chains = False
        
        # Extract total number of atoms
        atoms_match = re.search(r'(\d+) atoms', line)
        if atoms_match:
            total_atoms = int(atoms_match.group(1))
        
        # Extract total number of residues
        residues_match = re.search(r'(\d+) residues', line)
        if residues_match:
            total_residues = int(residues_match.group(1))

        # Check for missing atoms
        if "Missing heavy atoms: True" in line:
            missing_atoms = True

        # Check for multiple ligands
        ligand_match = re.search(rf"Multiple '{ligand_id}' ligand residues present: (True|False)", line)
        if ligand_match:
            multiple_ligands = ligand_match.group(1) == "True"

        # Check for multiple protein chains
        if "Multiple protein chains present: True" in line:
            multiple_protein_chains = True
    
    # Append the last entry
    if current_pdb:
        pdb_type = classify_pdb_type(total_atoms, total_residues, multiple_ligands, multiple_protein_chains)
        results.append((current_pdb, total_atoms, total_residues, missing_atoms, multiple_ligands, multiple_protein_chains, pdb_type))
    
    # Print results
    print(f"{'PDB File':<30} {'Total Atoms':<15} {'Total Residues':<15} {'Missing Atoms':<15} {'Multiple Ligands':<20} {'Multiple Protein Chains':<25} {'PDB Type':<15}")
    print("=" * 150)
    for pdb, atoms, residues, has_missing, has_multiple_ligands, has_multiple_chains, pdb_type in results:
        print(f"{pdb:<30} {str(atoms):<15} {str(residues):<15} {str(has_missing):<15} {str(has_multiple_ligands):<20} {str(has_multiple_chains):<25} {pdb_type:<15}")

def classify_pdb_type(total_atoms, total_residues, multiple_ligands, multiple_protein_chains):
    if total_atoms == 0 or total_residues == 0:
        return "not_available"
    if not multiple_ligands and not multiple_protein_chains:
        return "type_11"
    elif multiple_ligands and not multiple_protein_chains:
        return "type_n1"
    elif not multiple_ligands and multiple_protein_chains:
        return "type_1n"
    else:
        return "type_nn"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract missing atoms, multiple ligands, protein chain details, total atoms, total residues, and PDB type classification from a PDB log file.")
    parser.add_argument("input_file", type=str, help="Path to the log file containing PDB check outputs.")
    
    args = parser.parse_args()
    analyze_pdb_output(args.input_file)

