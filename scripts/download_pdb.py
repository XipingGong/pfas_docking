#!/usr/bin/env python3

import os
import requests
import argparse
import mdtraj as md
from Bio.PDB import MMCIFParser, PDBIO

def convert_mmcif_to_pdb(mmcif_filename, pdb_id, ligand_id=None, output_dir=".", keepcif=False):
    """
    Convert an mmCIF file to a PDB file using BioPython.

    Parameters:
        mmcif_filename (str): Path to the input mmCIF file.
        pdb_id (str): Desired PDB ID (used for naming the output file).
        ligand_id (str, optional): Ligand ID to be included in the output file name.
        output_dir (str, optional): Directory to save the output PDB file.
        keepcif (bool, optional): Whether to keep the mmCIF file after conversion.

    Returns:
        str or None: Path to the generated PDB file, or None if conversion failed.
    """
    filename_suffix = f"_{ligand_id}" if ligand_id else ""
    pdb_filename = os.path.join(output_dir, f"{pdb_id}{filename_suffix}.pdb")

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, mmcif_filename)

        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_filename)

        print(f"Converted {mmcif_filename} to {pdb_filename} successfully.")

        # Delete mmCIF file if not keeping it
        if not keepcif:
            os.remove(mmcif_filename)
            print(f"Deleted {mmcif_filename} as --keepcif is set to False.")

        return pdb_filename

    except Exception as e:
        print(f"Error converting {mmcif_filename} to PDB: {e}")

        # Delete mmCIF file even if conversion fails
        if not keepcif:
            os.remove(mmcif_filename)
            print(f"Deleted {mmcif_filename} due to conversion failure.")
        
        return None

def download_pdb(pdb_id, ligand_id=None, output_dir=".", overwrite=False, keepcif=False, cleanpdb=True):
    """
    Download a PDB file from RCSB and save it to the specified directory.
    If unavailable, attempt to download the mmCIF file and convert it to PDB.

    Parameters:
        pdb_id (str): PDB ID of the structure to download.
        ligand_id (str, optional): Ligand ID to include in the output file name.
        output_dir (str, optional): Directory to save the downloaded file.
        overwrite (bool, optional): Whether to overwrite an existing file.
        keepcif (bool, optional): Whether to keep the mmCIF file after conversion.
        cleanpdb (bool, optional): Whether to clean the downloaded PDB file using MDTraj.

    Returns:
        str or None: Path to the downloaded/converted PDB file, or None if failed.
    """
    os.makedirs(output_dir, exist_ok=True)
    filename_suffix = f"_{ligand_id}" if ligand_id else ""
    pdb_filename = os.path.join(output_dir, f"{pdb_id}{filename_suffix}.pdb")
    mmcif_filename = os.path.join(output_dir, f"{pdb_id}.cif")

    if os.path.exists(pdb_filename) and not overwrite:
        print(f"File {pdb_filename} already exists. Skipping download (overwrite=False).")
    else:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        mmcif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"

        response = requests.get(pdb_url)
        if response.status_code == 200:
            with open(pdb_filename, "w") as file:
                file.write(response.text)
            print(f"PDB file {pdb_filename} downloaded successfully.")
        else:
            print(f"PDB file not available for {pdb_id}, attempting to download mmCIF file.")
            response = requests.get(mmcif_url)
            if response.status_code == 200:
                with open(mmcif_filename, "w") as file:
                    file.write(response.text)
                print(f"mmCIF file {mmcif_filename} downloaded successfully.")
                pdb_filename = convert_mmcif_to_pdb(mmcif_filename, pdb_id, ligand_id, output_dir, keepcif)
            else:
                print(f"Failed to download both PDB and mmCIF files for ID {pdb_id}.")
                return None

    # Load PDB with MDTraj and print trajectory info
    try:
        traj = md.load(pdb_filename)
        print(traj)  # Prints trajectory info
    except Exception as e:
        print(f"Error loading {pdb_filename} with MDTraj: {e}")
        return pdb_filename

    # Optionally clean and save the PDB file
    if cleanpdb:
        cleaned_pdb_filename = pdb_filename  # Overwrite original file
        try:
            traj.save_pdb(cleaned_pdb_filename)
            print(f"Cleaned PDB file saved (without heading info): {cleaned_pdb_filename}")
        except Exception as e:
            print(f"Error cleaning and saving PDB file: {e}")

    return pdb_filename

def main():
    parser = argparse.ArgumentParser(description="Download a PDB or mmCIF file from RCSB, convert if necessary, and optionally clean with MDTraj.")
    parser.add_argument("pdb_id", type=str, help="PDB ID of the structure to download.")
    parser.add_argument("--ligand_id", type=str, default=None, help="Ligand ID to include in the output file name (default: None).")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the downloaded file (default: current directory).")
    parser.add_argument("--overwrite", type=lambda x: (str(x).lower() == 'true'), default=False, help="Whether to overwrite an existing file (default: False).")
    parser.add_argument("--keepcif", type=lambda x: (str(x).lower() == 'true'), default=False, help="Whether to keep the mmCIF file after conversion (default: False).")
    parser.add_argument("--cleanpdb", type=lambda x: (str(x).lower() == 'true'), default=True, help="Whether to clean the PDB file using MDTraj (default: True).")

    args = parser.parse_args()
    download_pdb(args.pdb_id, args.ligand_id, args.output_dir, args.overwrite, args.keepcif, args.cleanpdb)

if __name__ == "__main__":
    main()

