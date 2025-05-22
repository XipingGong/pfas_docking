import mdtraj as md
import os
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def extract_pdb(pdb_file, ligand_id):
    """
    Extracts protein residues closest to each ligand molecule with the given ligand ID,
    then selects the entire protein chain containing the closest residue.

    Args:
        pdb_file (str): Path to the input PDB file.
        ligand_id (str): The residue name of the ligand.
    """
    # Load the PDB file
    traj = md.load_pdb(pdb_file)
    topology = traj.topology
    print(f"Importing: {pdb_file}, {traj}")

    # Count residues
    num_protein = sum(1 for res in traj.topology.residues if res.is_protein)
    num_ligand = sum(1 for res in traj.topology.residues if res.name == ligand_id)
    print(f"The number of protein residues: {num_protein}")
    print(f"The number of ligand {ligand_id} residues: {num_ligand}")

    # Get PDB file name without extension
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_folder = os.path.dirname(pdb_file)  # Ensure extracted files go to the same folder

    # Select all ligand atoms
    ligand_atoms = topology.select(f"resname '{ligand_id}'")
    protein_atoms = topology.select("protein")

    # Get unique ligand residue indices
    ligand_residues = {topology.atom(idx).residue.index for idx in ligand_atoms}
    print(f"{ligand_id} Ligand residue index list: {ligand_residues}")

    # Iterate through each ligand instance
    for lig_res_index in ligand_residues:
        # Select atoms belonging to the specific ligand
        lig_atoms_list = topology.select(f"resid {lig_res_index}")
        if len(lig_atoms_list) == 0:
            print(f"Warning: No atoms found for ligand {lig_res_index}. Skipping...")
            continue
        lig_atoms = lig_atoms_list[0]  # Ensure we use the first atom

        # Compute pairwise distances between ligand atoms and protein atoms
        if len(protein_atoms) == 0:
            print(f"No protein atoms found for ligand {lig_res_index}. Skipping...")
            continue

        distances = md.compute_distances(traj, np.array([[lig_atoms, prot] for prot in protein_atoms if prot != lig_atoms]))

        # Identify the closest protein residue
        closest_prot_atom = protein_atoms[np.argmin(distances)]
        closest_residue = topology.atom(closest_prot_atom).residue
        closest_chain_id = closest_residue.chain.chain_id
        closest_chain_index = closest_residue.chain.index

        # Select all atoms in the same chain as the closest residue
        chain_atoms = topology.select(
            f"((protein and chainid {closest_chain_index}) or (resid {lig_res_index})) and element != H"
        )

        # Extract sub-trajectory for the ligand and protein chain
        sub_traj = traj.atom_slice(chain_atoms)

        # Save temporary structure
        temp_pdb = os.path.join(output_folder, f"temp_{base_name}_{ligand_id}_{lig_res_index}.pdb")
        sub_traj.save_pdb(temp_pdb)

        # Use PDBFixer to clean residue numbering and avoid MDTraj warnings
        fixer = PDBFixer(filename=temp_pdb)
        output_pdb = os.path.join(output_folder, f"{base_name}_{ligand_id}_{lig_res_index}_Protein_Chain{closest_chain_id}.pdb")
        with open(output_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        os.remove(temp_pdb)  # Clean up temporary file
        print(f"Saved: {output_pdb}, {sub_traj}")

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3 or sys.argv[1] == '-h':
        print("Usage: python extract_pdb.py <pdb_file> <ligand_id>")
        print("Description:")
        print("  This script extracts the closest protein chain to a ligand in a given PDB file.")
        print("  It identifies ligand instances, finds the nearest protein residue, and saves")
        print("  the protein chain that contains this residue along with the ligand.")
        print("Arguments:")
        print("  <pdb_file>   Path to the input PDB file.")
        print("  <ligand_id>  The residue name of the ligand.")
        print("Example:")
        print("  python extract_pdb.py my_structure.pdb LIG")
        sys.exit(1)

    pdb_filename = sys.argv[1]
    ligand_name = sys.argv[2]
    extract_pdb(pdb_filename, ligand_name)

