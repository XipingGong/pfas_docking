#!/usr/bin/env python3

import mdtraj as md
import argparse

def detect_ligand_id(traj):
    """
    Detect the ligand ID by finding non-protein residues.

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

    return list(ligand_resnames)[0]

def save_selected_atoms_from_pdb(pdb_file, selected_indices, output_file):
    """
    Save selected atoms to a PDB file without altering the original format.

    Args:
        pdb_file (str): Path to the original PDB file.
        selected_indices (set[int]): Set of atom indices (zero-based) to save.
        output_file (str): Output path for the selected atoms.
    """
    selected_indices = set(selected_indices)
    current_atom_index = 0  # Tracks ATOM/HETATM lines

    with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                if current_atom_index in selected_indices:
                    f_out.write(line)
                current_atom_index += 1
           #elif line.startswith('TER') or line.startswith('END'):
           #    f_out.write(line)

def get_inputs_for_vina(pdb_file, pocket_size_ang=20,
                        protein_output="x_receptor.pdb",
                        ligand_output="x_ligand.pdb",
                        pocket_params_output="x_pocket_params.txt"):
    """
    Extract protein and ligand from a PDB file and compute docking pocket parameters.

    Args:
        pdb_file (str): Path to input PDB file.
        pocket_size_ang (float): Size of the docking pocket in Ångstroms.
        protein_output (str): Output path for protein structure.
        ligand_output (str): Output path for ligand structure.
        pocket_params_output (str): Output path for docking pocket parameters.
    """
    try:
        print(f"📂 Loading PDB file: {pdb_file}")
        traj = md.load(pdb_file)
        topology = traj.topology

        ligand_id = detect_ligand_id(traj)
        print(f"🔍 Detected ligand ID: {ligand_id}")

        protein_selection = topology.select("protein")
        if protein_selection.size == 0:
            raise ValueError("❌ Error: No protein found in the structure.")

        ligand_selection = topology.select(f"resname '{ligand_id}'")
        if ligand_selection.size == 0:
            raise ValueError(f"❌ Error: No ligand with residue name '{ligand_id}' found.")

        # Compute docking pocket center
        pocket_center = md.compute_center_of_mass(traj.atom_slice(ligand_selection))[0]
        pocket_center_ang = pocket_center * 10  # nm -> Å

        center_x, center_y, center_z = pocket_center_ang
        size_x = size_y = size_z = pocket_size_ang

        # Save pocket parameters
        with open(pocket_params_output, "w") as f:
            f.write(f"center_x = {center_x:.3f}\n")
            f.write(f"center_y = {center_y:.3f}\n")
            f.write(f"center_z = {center_z:.3f}\n")
            f.write(f"size_x = {size_x:.3f}\n")
            f.write(f"size_y = {size_y:.3f}\n")
            f.write(f"size_z = {size_z:.3f}\n")

        print(f"📌 Docking pocket parameters saved to: {pocket_params_output}")

        # Save protein and ligand separately
        save_selected_atoms_from_pdb(pdb_file, protein_selection, protein_output)
        save_selected_atoms_from_pdb(pdb_file, ligand_selection, ligand_output)

        print(f"✅ Successfully extracted:\n"
              f"   - Protein: {protein_output}\n"
              f"   - Ligand: {ligand_output}\n"
              f"   - Docking Box Params: {pocket_params_output}")

    except Exception as e:
        print(f"❌ Error: {str(e)}")

def main():
    """Main function to handle argument parsing."""
    parser = argparse.ArgumentParser(
        description="Extract protein and ligand from a PDB file and compute docking pocket parameters for AutoDock Vina."
    )
    parser.add_argument("pdb_file", type=str, help="Path to input PDB file.")
    parser.add_argument("--pocket_size_ang", type=float, default=20,
                        help="Size of the docking pocket in Ångstroms (default: 20).")
    parser.add_argument("--protein_output", type=str, default="x_receptor.pdb",
                        help="Output path for protein structure (default: x_receptor.pdb).")
    parser.add_argument("--ligand_output", type=str, default="x_ligand.pdb",
                        help="Output path for ligand structure (default: x_ligand.pdb).")
    parser.add_argument("--pocket_params_output", type=str, default="x_pocket_params.txt",
                        help="Output path for docking pocket parameters (default: x_pocket_params.txt).")

    args = parser.parse_args()

    get_inputs_for_vina(
        args.pdb_file,
        args.pocket_size_ang,
        args.protein_output,
        args.ligand_output,
        args.pocket_params_output
    )

if __name__ == "__main__":
    main()

