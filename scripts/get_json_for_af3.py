import argparse
import json
import mdtraj as md
import os

def extract_protein_sequence_and_ligand(pdb_filename):
    """
    Extracts the protein sequence and ligand ID from a PDB file using MDTraj.

    Parameters:
    - pdb_filename (str): Path to the input PDB file.

    Returns:
    - dict: JSON-compatible dictionary with protein sequence and ligand info.
    """
    # Load the PDB file using MDTraj
    traj = md.load(pdb_filename)

    # Extract protein sequences using MDTraj's built-in function
    fasta_sequences = traj.topology.to_fasta() if traj.topology.n_residues > 0 else []

    # Clean up empty sequences (some chains may be empty)
    fasta_sequences = [seq for seq in fasta_sequences if seq.strip()]

    # Extract chain IDs from the PDB file (not using `index`)
    chain_ids = []
    for chain in traj.topology.chains:
        first_residue = next((res for res in chain.residues if res.is_protein), None)
        if first_residue:
            chain_ids.append(first_residue.chain.chain_id)

    # Ensure protein sequences and chain IDs match
    if fasta_sequences and len(fasta_sequences) != len(chain_ids):
        raise ValueError("Mismatch between protein sequences and chain IDs. Check the PDB file.")

    # Create JSON entries for protein sequences
    protein_sequences = [
        {"protein": {"id": chain_ids[i], "sequence": sequence}}
        for i, sequence in enumerate(fasta_sequences) if sequence
    ]

    # Extract ligand information
    ligands = []
    for residue in traj.topology.residues:
        if not residue.is_protein:
            ligand_info = {"id": residue.chain.chain_id, "ccdCodes": [residue.name]}
            ligands.append(ligand_info)

    # Ensure at least one valid sequence or ligand exists
    if not protein_sequences and not ligands:
        raise ValueError("No valid protein sequence or ligand found in the PDB file.")

    # Prepare JSON output
    json_data = {
        "name": os.path.splitext(os.path.basename(pdb_filename))[0],  # Use filename (without .pdb) as name
        "sequences": protein_sequences + ([{"ligand": ligands[0]}] if ligands else []),  
        "modelSeeds": [1],
        "bondedAtomPairs": [],
        "dialect": "alphafold3",
        "version": 2
    }

    return json_data

def save_json(data, output_filename):
    """Save extracted PDB info to a JSON file."""
    with open(output_filename, "w") as json_file:
        json.dump(data, json_file, indent=2)

def main():
    """Main function to handle argument parsing and execution."""
    parser = argparse.ArgumentParser(description="Extract protein sequences and ligand information from a PDB file and save as JSON.")
    
    parser.add_argument("input_pdb", help="Path to the input PDB file (1 protein, chain A & 1 ligand, chain B, and ligand_id is the ccdCode.)")
    parser.add_argument("-o", "--output", help="Path to save the output JSON file (default: <input_pdb>.json).")

    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_json = args.output if args.output else input_pdb.replace(".pdb", ".json")

    try:
        extracted_data = extract_protein_sequence_and_ligand(input_pdb)
        save_json(extracted_data, output_json)
        print(f"✅ JSON file saved: {output_json}")
    except Exception as e:
        print(f"❌ Error: {e}")
        exit(1)

if __name__ == "__main__":
    main()

