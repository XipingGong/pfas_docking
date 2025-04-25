#!/usr/bin/env python3

import os
import requests
import argparse
import mdtraj as md

def download_file(url, save_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "w") as f:
            f.write(response.text)
        return True
    return False

def download_and_convert(pdb_id, ligand_id=None, output_dir=".", overwrite=False):
    """
    Download a structure file (PDB or CIF), load with MDTraj, and save as cleaned PDB.
    """
    os.makedirs(output_dir, exist_ok=True)
    suffix = f"_{ligand_id}" if ligand_id else ""
    output_pdb_path = os.path.join(output_dir, f"{pdb_id}{suffix}.pdb")

    if os.path.exists(output_pdb_path) and not overwrite:
        print(f"✅ PDB already exists: {output_pdb_path} (use --overwrite to re-download)")
        return output_pdb_path

    # Try downloading PDB first
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    tmp_file = os.path.join(output_dir, f"{pdb_id}_raw.pdb")

    if not download_file(pdb_url, tmp_file):
        print(f"⚠️ PDB not available for {pdb_id}, trying mmCIF...")
        cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        tmp_file = os.path.join(output_dir, f"{pdb_id}.cif")

        if not download_file(cif_url, tmp_file):
            print(f"❌ Failed to download both PDB and mmCIF for {pdb_id}.")
            return None
        else:
            print(f"📥 Downloaded mmCIF: {tmp_file}")
    else:
        print(f"📥 Downloaded PDB: {tmp_file}")

    # Try loading and saving as clean PDB
    try:
        traj = md.load(tmp_file)
        traj.save_pdb(output_pdb_path)
        print(f"✅ Saved cleaned PDB: {output_pdb_path}")
        return output_pdb_path
    except Exception as e:
        print(f"❌ MDTraj failed to process {tmp_file}: {e}")
        return None
    finally:
        # Clean up temporary raw files
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

def main():
    parser = argparse.ArgumentParser(description="Download a structure and convert to cleaned PDB using MDTraj.")
    parser.add_argument("pdb_id", type=str, help="PDB ID to download.")
    parser.add_argument("--ligand_id", type=str, default=None, help="Ligand ID for naming the output.")
    parser.add_argument("--output_dir", type=str, default=".", help="Directory to save the final PDB.")
    parser.add_argument("--overwrite", type=lambda x: str(x).lower() == 'true', default=False, help="Overwrite existing PDB if it exists.")

    args = parser.parse_args()
    download_and_convert(args.pdb_id, args.ligand_id, args.output_dir, args.overwrite)

if __name__ == "__main__":
    main()

