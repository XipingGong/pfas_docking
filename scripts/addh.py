#!/usr/bin/env python3
"""
addh.py: Update heavy and hydrogen atom coordinates in a reference PDB
using an input heavy-only PDB and a reference MOL2 file.
If the MOL2 does not contain heavy-hydrogen bonds, fall back to distance-based guessing.

Usage:
  python addh.py input_heavy.pdb --ref_pdb ref_with_H.pdb --ref_mol2 ref_with_H.mol2 --output_pdb output_with_H.pdb
"""

import argparse
import numpy as np
import sys


def is_atom_line(line):
    return line.startswith(('ATOM  ', 'HETATM'))


def parse_input_pdb(lines):
    heavy_map = {}
    first_res_seq = None
    for line in lines:
        if not is_atom_line(line):
            continue
        atom_name = line[12:16].strip()
        if atom_name.startswith('H'):
            continue
        res_name = line[17:20].strip()
        chain = line[21]
        res_seq = int(line[22:26].strip())
        i_code = line[26]

        if first_res_seq is None:
            first_res_seq = res_seq

        rel_res_seq = str(res_seq - first_res_seq)
        key = (chain, rel_res_seq, i_code, res_name, atom_name)

        raw = line[30:54]
        x = float(raw[0:8])
        y = float(raw[8:16])
        z = float(raw[16:24])
        heavy_map[key] = {
            'raw': raw,
            'xyz': np.array([x, y, z])
        }
    return heavy_map


def parse_ref_pdb(lines):
    entries = []
    first_res_seq = None
    for line in lines:
        if is_atom_line(line):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_seq = int(line[22:26].strip())
            i_code = line[26]

            if first_res_seq is None:
                first_res_seq = res_seq

            rel_res_seq = str(res_seq - first_res_seq)
            key = (chain, rel_res_seq, i_code, res_name, atom_name)

            entries.append({
                'key': key,
                'orig_line': line,
                'new_raw': line[30:54],
                'new_xyz': None
            })
    return entries, lines


def parse_mol2_atoms_and_bonds(lines):
    atoms = {}
    bonds = []
    atom_sec = bond_sec = False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            atom_sec, bond_sec = True, False
            continue
        if line.startswith("@<TRIPOS>BOND"):
            atom_sec, bond_sec = False, True
            continue
        if line.startswith("@<TRIPOS>"):
            atom_sec = bond_sec = False
        if atom_sec:
            parts = line.split()
            aid = int(parts[0])
            coords_ref = np.array(list(map(float, parts[2:5])))
            atoms[aid] = {
                'coords_ref': coords_ref,
                'type': parts[5]
            }
        if bond_sec:
            parts = line.split()
            bonds.append((int(parts[1]), int(parts[2])))
    return atoms, bonds


def mol2_has_complete_hydrogen_bonds(mol2_atoms, bonds):
    bonded_hydrogens = set()
    heavy_atom_ids = {aid for aid, atom in mol2_atoms.items() if not atom['type'].startswith('H')}

    for i, j in bonds:
        if i in mol2_atoms and j in mol2_atoms:
            if mol2_atoms[i]['type'].startswith('H') and j in heavy_atom_ids:
                bonded_hydrogens.add(i)
            if mol2_atoms[j]['type'].startswith('H') and i in heavy_atom_ids:
                bonded_hydrogens.add(j)

    all_hydrogens = {aid for aid, atom in mol2_atoms.items() if atom['type'].startswith('H')}
    missing_hydrogens = all_hydrogens - bonded_hydrogens

    if missing_hydrogens:
        print(f"Warning: {len(missing_hydrogens)} hydrogens have no bonded heavy atom in MOL2.")
        return False
    return True


def update_heavy(entries, heavy_map):
    for entry in entries:
        key = entry['key']
        atom_name = key[-1]
        if key in heavy_map:
            info = heavy_map[key]
            entry['new_raw'] = info['raw']
            entry['new_xyz'] = info['xyz']
        elif not atom_name.startswith('H'):
            sys.exit(f"ERROR: Cannot find the heavy atom: {key} in input_pdb")
    return entries


def project_hydrogens_mol2(entries, mol2_atoms, bonds):
    for aid, atom in mol2_atoms.items():
        if atom['type'].startswith('H'):
            idx = aid - 1
            if idx >= len(entries):
                continue
            entry = entries[idx]
            partner = None
            for i, j in bonds:
                if aid == i and not mol2_atoms[j]['type'].startswith('H'):
                    partner = j
                    break
                if aid == j and not mol2_atoms[i]['type'].startswith('H'):
                    partner = i
                    break
            if partner is None:
                sys.exit(f"ERROR: No heavy partner for H id {aid}")

            heavy_entry = entries[partner - 1]
            if heavy_entry['new_xyz'] is None:
                sys.exit(f"ERROR: Missing heavy atom coordinates for bonded H id {aid}")

            vec_ref = atom['coords_ref'] - mol2_atoms[partner]['coords_ref']
            H_new = heavy_entry['new_xyz'] + vec_ref
            x, y, z = H_new
            entry['new_raw'] = f"{x:8.3f}{y:8.3f}{z:8.3f}"
    return entries


def project_hydrogens_by_distance(entries, cutoff=2.0):
    heavy_entries = [e for e in entries if not e['key'][-1].startswith('H')]
    hydrogen_entries = [e for e in entries if e['key'][-1].startswith('H')]

    heavy_positions = np.array([
        np.array([
            float(e['orig_line'][30:38]),
            float(e['orig_line'][38:46]),
            float(e['orig_line'][46:54])
        ]) for e in heavy_entries
    ])
    heavy_positions_new = np.array([e['new_xyz'] for e in heavy_entries])

    for h_entry in hydrogen_entries:
        h_pos_old = np.array([
            float(h_entry['orig_line'][30:38]),
            float(h_entry['orig_line'][38:46]),
            float(h_entry['orig_line'][46:54])
        ])

        distances = np.linalg.norm(heavy_positions - h_pos_old, axis=1)
        nearest_idx = np.argmin(distances)
        nearest_dist = distances[nearest_idx]

        if nearest_dist > cutoff:
            print(f"Warning: Hydrogen {h_entry['key']} is far from heavy atoms ({nearest_dist:.2f} Å). Keeping original position.")
            continue

        delta = heavy_positions_new[nearest_idx] - heavy_positions[nearest_idx]
        h_pos_new = h_pos_old + delta
        x, y, z = h_pos_new
        h_entry['new_raw'] = f"{x:8.3f}{y:8.3f}{z:8.3f}"

    return entries


def write_output(ref_lines, entries, out_path):
    atom_idx = 0
    with open(out_path, 'w') as fout:
        for line in ref_lines:
            if is_atom_line(line):
                entry = entries[atom_idx]
                new_line = line[:30] + entry['new_raw'] + line[54:]
                fout.write(new_line)
                atom_idx += 1
            else:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_pdb', help='Heavy-only PDB (input)')
    parser.add_argument('--ref_pdb', required=True, help='Reference PDB (with H)')
    parser.add_argument('--ref_mol2', required=True, help='Reference MOL2 (with H, bonds)')
    parser.add_argument('-o', '--output_pdb', required=True, help='Output PDB with updated coordinates')
    args = parser.parse_args()

    with open(args.input_pdb) as f:
        input_lines = f.readlines()
    with open(args.ref_pdb) as f:
        ref_lines = f.readlines()
    with open(args.ref_mol2) as f:
        mol2_lines = f.readlines()

    heavy_map = parse_input_pdb(input_lines)
    entries, ref_all_lines = parse_ref_pdb(ref_lines)
    mol2_atoms, bonds = parse_mol2_atoms_and_bonds(mol2_lines)

    entries = update_heavy(entries, heavy_map)

    if mol2_has_complete_hydrogen_bonds(mol2_atoms, bonds):
        print("✅ Using MOL2 bond-based hydrogen projection...")
        entries = project_hydrogens_mol2(entries, mol2_atoms, bonds)
    else:
        print("⚠️  Falling back to distance-based hydrogen projection (cutoff = 2.0 Å)...")
        entries = project_hydrogens_by_distance(entries, cutoff=2.0)

    write_output(ref_all_lines, entries, args.output_pdb)


if __name__ == '__main__':
    main()

