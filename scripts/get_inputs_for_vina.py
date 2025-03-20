#!/usr/bin/env python3

import mdtraj as md
import os
import argparse

def detect_ligand_id(traj):
    """
    Detects the ligand ID from the PDB file by finding non-protein residues.

    Args:
        traj (md.Trajectory): MDTraj trajectory object.

    Returns:
        str: Detected ligand residue name.
    """
    ligand_resnames = {res.name for res in traj.topology.residues if not res.is_protein}
    
    if not ligand_resnames:
        raise ValueError("❌ Error: No ligand found in the PDB file.")
    
    if len(ligand_resnames) > 1:
        print(f"⚠ Warning: Multiple ligands detected: {ligand_resnames}. Using the first one.")

    return list(ligand_resnames)[0]  # Return the first detected ligand ID


def get_inputs_for_vina(pdb_file, pocket_size_ang=20, protein_output="x_receptor.pdb", 
                        ligand_output="x_ligand.pdb", box_params_output="x_box_params.txt"):
    """
    Extracts the protein and ligand structures from a PDB file and computes docking box parameters.

    Args:
        pdb_file (str): Path to the input PDB file.
        pocket_size_ang (float): Size of the docking box in Ångstroms.
        protein_output (str): Output path for protein structure.
        ligand_output (str): Output path for ligand structure.
        box_params_output (str): Output path for box parameters file.
    """
    try:
        print(f"📂 Loading PDB file: {pdb_file}")
        traj = md.load(pdb_file)
        topology = traj.topology

        # Automatically detect ligand ID
        ligand_id = detect_ligand_id(traj)
        print(f"🔍 Detected ligand ID: {ligand_id}")

        # Select protein atoms (standard protein residues)
        protein_selection = topology.select("protein")

        if protein_selection.size == 0:
            raise ValueError("❌ Error: No protein found in the structure.")

        # Select ligand atoms by residue name
        ligand_selection = topology.select(f"resname '{ligand_id}'")

        if ligand_selection.size == 0:
            raise ValueError(f"❌ Error: No ligand with residue name '{ligand_id}' found.")

        # Compute the docking box center and size
        pocket_center = md.compute_center_of_mass(traj.atom_slice(ligand_selection))[0]
        pocket_center_ang = pocket_center * 10  # Convert from nm to Å

        center_x, center_y, center_z = pocket_center_ang
        size_x = size_y = size_z = pocket_size_ang

        # Save box parameters to a file
        with open(box_params_output, "w") as f:
            f.write(f"center_x = {center_x:.3f}\n")
            f.write(f"center_y = {center_y:.3f}\n")
            f.write(f"center_z = {center_z:.3f}\n")
            f.write(f"size_x = {size_x:.3f}\n")
            f.write(f"size_y = {size_y:.3f}\n")
            f.write(f"size_z = {size_z:.3f}\n")

        print(f"📌 Docking box parameters saved to: {box_params_output}")

        # Extract and save protein and ligand structures
        protein_traj = traj.atom_slice(protein_selection)
        ligand_traj = traj.atom_slice(ligand_selection)

        protein_traj.save(protein_output)
        ligand_traj.save(ligand_output)

        print(f"✅ Successfully extracted:\n"
              f"   - Protein: {protein_output}\n"
              f"   - Ligand: {ligand_output}\n"
              f"   - Docking Box Params: {box_params_output}")

    except Exception as e:
        print(f"❌ Error: {str(e)}")


def main():
    """Main function to handle argument parsing."""
    parser = argparse.ArgumentParser(
        description="Extract protein and ligand from a PDB file and compute docking box parameters for AutoDock Vina."
    )
    parser.add_argument("pdb_file", type=str, help="Path to input PDB file.")
    parser.add_argument("--pocket_size_ang", type=float, default=20, help="Size of the docking box in Ångstroms (default: 20).")
    parser.add_argument("--protein_output", type=str, default="x_receptor.pdb", help="Output path for protein structure (default: x_receptor.pdb).")
    parser.add_argument("--ligand_output", type=str, default="x_ligand.pdb", help="Output path for ligand structure (default: x_ligand.pdb).")
    parser.add_argument("--box_params_output", type=str, default="x_box_params.txt", help="Output path for docking box parameters (default: x_box_params.txt).")

    args = parser.parse_args()

    get_inputs_for_vina(
        args.pdb_file, args.pocket_size_ang, args.protein_output, args.ligand_output, args.box_params_output
    )


if __name__ == "__main__":
    main()

