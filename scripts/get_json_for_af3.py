import argparse
import json
import mdtraj as md
import os

def extract_and_save_json(pdb_filename, output_filename):
    """
    Extract protein sequence and ligand info from a PDB file and save as JSON.

    Parameters:
    - pdb_filename (str): Path to the input PDB file.
    - output_filename (str): Path to the output JSON file.
    """
    traj = md.load(pdb_filename)

    fasta_sequences = traj.topology.to_fasta() if traj.topology.n_residues > 0 else []
    fasta_sequences = [seq for seq in fasta_sequences if seq.strip()]

    chain_ids = []
    for chain in traj.topology.chains:
        first_residue = next((res for res in chain.residues if res.is_protein), None)
        if first_residue:
            chain_ids.append(first_residue.chain.chain_id)

    if fasta_sequences and len(fasta_sequences) != len(chain_ids):
        raise ValueError("Mismatch between protein sequences and chain IDs. Check the PDB file.")

    protein_sequences = [
        {"protein": {"id": chain_ids[i], "sequence": sequence}}
        for i, sequence in enumerate(fasta_sequences) if sequence
    ]

    ligands = []
    for residue in traj.topology.residues:
        if not residue.is_protein:
            ligand_info = {"id": residue.chain.chain_id, "ccdCodes": [residue.name]}
            ligands.append(ligand_info)

    if not protein_sequences and not ligands:
        raise ValueError("No valid protein sequence or ligand found in the PDB file.")

    json_data = {
        "name": os.path.splitext(os.path.basename(output_filename))[0],
        "sequences": protein_sequences + ([{"ligand": ligands[0]}] if ligands else []),
        "modelSeeds": [1],
        "bondedAtomPairs": [],
        "dialect": "alphafold3",
        "version": 2
    }

    with open(output_filename, "w") as json_file:
        json.dump(json_data, json_file, indent=2)

def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences and ligand info from a PDB file and save as JSON.")
    parser.add_argument("input_pdb", help="Path to the input PDB file (1 protein, chain A & 1 ligand, chain B, and ligand_id is the ccdCode.)")
    parser.add_argument("-o", "--output", help="Path to save the output JSON file (default: <input_pdb>.json).")

    args = parser.parse_args()
    input_pdb = args.input_pdb
    output_json = args.output if args.output else input_pdb.replace(".pdb", ".json")

    try:
        extract_and_save_json(input_pdb, output_json)
        print(f"✅ JSON file saved: {output_json}")
    except Exception as e:
        print(f"❌ Error: {e}")
        exit(1)

if __name__ == "__main__":
    main()

