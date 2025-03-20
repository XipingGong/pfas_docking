import numpy as np
import mdtraj as md
import sys
import os
import argparse
import datetime

def log_message(message):
    """Logs a timestamped message."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def analyze_ligand_interactions(traj, cutoff_distance, ligand_resid, traj_name = ''):
    """
    Analyze residues near a ligand within a cutoff distance and categorize them by type.

    Parameters:
    - traj: mdtraj.Trajectory
        The loaded molecular dynamics trajectory.
    - cutoff_distance: float
        The cutoff distance (in nm) to find nearby residues.
    - ligand_resid: int
        The MDTraj-supported residue index of the ligand (0-based).
    - traj_name: string
        The name of traj.

    Returns:
    - None
        Prints categorized residues and outputs a selection string for VMD visualization.
    """
    # Extract the chain number and PDB residue ID (resSeq) from the residue index
    residue = traj.topology.residue(ligand_resid)
    ligand_chain_id = residue.chain.chain_id
    ligand_resSeq = residue.resSeq
    ligand_name = residue.name
    if not ligand_name:
        print(f"Ligand with residue ID {ligand_resid} not found!")
        return

    print(f"# Identifying the residues within {cutoff_distance} nm of a ligand (chain_id: {ligand_chain_id}, resSeq: {ligand_resSeq})")
    print("# ------------------------------------------------------------------------------------")

    # Define residue groups based on properties
    non_polar = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
    polar = {"SER", "THR", "CYS", "TYR", "ASN", "GLN"}
    positively_charged = {"LYS", "ARG", "HIS"}
    negatively_charged = {"ASP", "GLU"}
    categories = {"Non-Polar": [], "Polar": [], "Positively Charged": [], "Negatively Charged": []}

    # Find atoms within the cutoff distance of the ligand
    nearby_atoms = md.compute_neighbors(
        traj,
        cutoff_distance,
        query_indices=traj.topology.select(f"resid {ligand_resid} and (not element H)"),
        haystack_indices=traj.topology.select("protein and (not element H)"),
        periodic=True
    )
    # print(f"nearby_atoms = {nearby_atoms}")

    # Get the residue IDs and names of nearby atoms
    residue_info = set(
        (traj.topology.atom(idx).residue.chain.chain_id,
         traj.topology.atom(idx).residue.index,
         traj.topology.atom(idx).residue.resSeq,
         traj.topology.atom(idx).residue.name)
        for idx in nearby_atoms[0]
    )

    # Categorize residues
    for chain_id, resid, resSeq, resName in residue_info:
        if resName in non_polar:
            categories["Non-Polar"].append((chain_id, resid, resSeq, resName))
        elif resName in polar:
            categories["Polar"].append((chain_id, resid, resSeq, resName))
        elif resName in positively_charged:
            categories["Positively Charged"].append((chain_id, resid, resSeq, resName))
        elif resName in negatively_charged:
            categories["Negatively Charged"].append((chain_id, resid, resSeq, resName))

    # Print categorized residues and their distances
    for group, residues in categories.items():
        print(f"{group} Residues:")
        if residues:
            for _, resid, _, _ in residues:
                # Calculate the minimum distance between ligand and residue
                residue1_atoms = traj.topology.select(f"resid {ligand_resid} and (not element H)")
                residue2_atoms = traj.topology.select(f"resid {resid} and (not element H)")
                pairs = np.array([(a1, a2) for a1 in residue1_atoms for a2 in residue2_atoms])

                distances = md.compute_distances(traj, pairs)
                min_distance = distances.min()
                min_index = distances.argmin()

                chain_id_1 = traj.topology.atom(pairs[min_index][0]).residue.chain.chain_id
                resid_1    = traj.topology.atom(pairs[min_index][0]).residue.index
                resSeq_1   = traj.topology.atom(pairs[min_index][0]).residue.resSeq
                resname_1  = traj.topology.atom(pairs[min_index][0]).residue.name
                atomtype_1 = traj.topology.atom(pairs[min_index][0]).name

                chain_id_2 = traj.topology.atom(pairs[min_index][1]).residue.chain.chain_id
                resid_2    = traj.topology.atom(pairs[min_index][1]).residue.index
                resSeq_2   = traj.topology.atom(pairs[min_index][1]).residue.resSeq
                resname_2  = traj.topology.atom(pairs[min_index][1]).residue.name
                atomtype_2 = traj.topology.atom(pairs[min_index][1]).name


                # Print the residue and distance information
                print(
                    f"  {traj_name} - Ligand (chain_id: {chain_id_1}, resSeq: {resSeq_1}, resName: {resname_1}), "
                    f"Residue (chain_id: {chain_id_2}, resSeq: {resSeq_2}, resName: {resname_2}), "
                    f"Min-Dist ({resname_1}-{atomtype_1},{resname_2}-{atomtype_2}): {min_distance:.2f} nm"
                    f"{' ***' if min_distance < 0.3 else ''}"
                )

        else:
            print("  N/A")

    # Print the VMD selection string
    ligand_string = f"(chain {ligand_chain_id} and resname {ligand_name} and resid {ligand_resSeq})"
    residues_string = " or ".join([f"(chain {chain_id} and resid {resSeq})" for chain_id, resid, resSeq, resName in residue_info])
    print(f"\nFor VMD visualization: {ligand_string} or {residues_string} \n")

def main():
    """Main function to process a PDB file and analyze ligand interactions."""
    parser = argparse.ArgumentParser(description="Analyze ligand interactions in a PDB file using MDTraj.")
    parser.add_argument("pdb_file", type=str, help="Path to the PDB file.")
    parser.add_argument("ligand_name", type=str, help="Ligand name to analyze.")
    parser.add_argument("--cutoff_distance", type=float, default=0.5, help="Cutoff distance in nm for identifying nearby residues (default: 0.5 nm).")
    args = parser.parse_args()

    pdb_file = args.pdb_file
    ligand_name = args.ligand_name
    cutoff_distance = args.cutoff_distance

    log_message(f"Processing PDB file: {pdb_file}, Ligand: {ligand_name}")

    # Load the trajectory
    log_message(f"Loading trajectory from {pdb_file}...")
    try:
        traj = md.load(pdb_file)
        #traj = traj[0]  # only using first model
        log_message(f"Trajectory loaded with {traj.n_atoms} atoms and {traj.n_residues} residues.")
    except Exception as e:
        log_message(f"[ERROR] Failed to load PDB file: {e}")
        sys.exit(1)

    log_message(f"Using cutoff distance: {cutoff_distance} nm")

    residue_ids = [res.index for res in traj.topology.residues if res.name == ligand_name]
    if not residue_ids:
        log_message(f"[WARNING] No residues found for ligand {ligand_name}. Check if the ligand name is correct.")
        sys.exit(1)

    for ligand_resid in residue_ids:
        log_message(f"Analyzing interactions for ligand residue ID: {ligand_resid}")
        try:
            analyze_ligand_interactions(traj, cutoff_distance, ligand_resid, traj_name=os.path.basename(pdb_file).split('.')[0])
        except Exception as e:
            log_message(f"[ERROR] Failed to analyze interactions for residue {ligand_resid}: {e}")

    log_message("Processing complete.")

if __name__ == "__main__":
    main()
