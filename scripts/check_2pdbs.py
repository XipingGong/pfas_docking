#!/usr/bin/env python3

import argparse
import numpy as np
import mdtraj as md

def parse_pdb_for_identity(pdb_path):
    records = []
    with open(pdb_path, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                identity = line[:30] + line[54:]
                records.append(identity.strip())
    return records

def calculate_rmsd(pdb1, pdb2):
    traj1 = md.load(pdb1)
    traj2 = md.load(pdb2)

    # 🧠 Select non-protein heavy atoms
    heavy_indices = traj1.topology.select("not protein and not element H")
    ref_slice = traj1.atom_slice(heavy_indices)
    target_slice = traj2.atom_slice(heavy_indices)

    # 🚫 Direct (unaligned) RMSD
    diff = target_slice.xyz - ref_slice.xyz[0]
    rmsd_direct = np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=1))
    rmsd_direct_min = rmsd_direct.min()

    # ✅ Aligned RMSD using MDTraj
    rmsd_aligned = md.rmsd(target_slice, ref_slice)
    rmsd_aligned_min = rmsd_aligned.min()

    return rmsd_aligned_min, rmsd_direct_min

def main():
    parser = argparse.ArgumentParser(
        description="Compare two PDBs: check identity (excluding coordinates) and RMSD of heavy ligand atoms."
    )
    parser.add_argument("pdb1", help="First PDB file")
    parser.add_argument("pdb2", help="Second PDB file")
    args = parser.parse_args()

    try:
        id1 = parse_pdb_for_identity(args.pdb1)
        id2 = parse_pdb_for_identity(args.pdb2)

        if id1 == id2:
            rmsd_aligned, rmsd_direct = calculate_rmsd(args.pdb1, args.pdb2)
            print(f"✅ {args.pdb1} - {args.pdb2}: Same except coordinates. "
                  f"Ligand_RMSD_MDTraj = {rmsd_aligned:.3f} nm, Ligand_RMSD_Direct = {rmsd_direct:.3f} nm")
        else:
            print(f"❌ {args.pdb1} - {args.pdb2}: Differ beyond coordinates.")
            print("Differences (ignoring coordinates):")
            print("=" * 60)
            for i, (line1, line2) in enumerate(zip(id1, id2), 1):
                if line1 != line2:
                    print(f"Line {i}:")
                    print(f"  {args.pdb1}: {line1}")
                    print(f"  {args.pdb2}: {line2}")
                    print("-" * 60)
                    if i >= 2:
                        print("   ...")
                        break
            if len(id1) != len(id2):
                print(f"\n⚠️ Warning: Files have different number of ATOM/HETATM lines.")
                print(f"{args.pdb1}: {len(id1)} lines")
                print(f"{args.pdb2}: {len(id2)} lines")

    except FileNotFoundError as e:
        print(f"❌ File not found: {e.filename}")
    except Exception as e:
        print(f"❌ An error occurred: {e}")

if __name__ == "__main__":
    main()

