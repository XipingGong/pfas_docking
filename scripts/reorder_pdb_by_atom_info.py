#!/usr/bin/env python3
import argparse

def build_residue_renumbering_map(pdb_lines):
    """Assign a relative index to each unique residue based on its appearance."""
    seen = {}
    res_map = {}
    index = 0
    for line in pdb_lines:
        if line.startswith(("ATOM", "HETATM")):
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            key = (res_name, chain_id, res_seq)
            if key not in seen:
                seen[key] = index
                index += 1
            res_map[key] = seen[key]
    return res_map

def extract_atom_key_from_line(line, res_map):
    """Generate an atom matching key using relative residue index."""
    atom_name = line[12:16].strip()
    res_name = line[17:20].strip()
    chain_id = line[21].strip()
    res_seq = line[22:26].strip()
    rel_res_index = res_map.get((res_name, chain_id, res_seq))
    if rel_res_index is None:
        raise ValueError(f"Residue not found in map: {(res_name, chain_id, res_seq)}")
    return (atom_name, res_name, chain_id, rel_res_index)

def parse_pdb_atoms(filepath):
    """Extract ATOM/HETATM lines from a PDB file."""
    with open(filepath) as f:
        return [line for line in f if line.startswith(('ATOM', 'HETATM'))]

def build_coord_line_map(pdb_lines, res_map):
    """Build a map from atom key → coordinate string (preserving exact formatting)."""
    coord_map = {}
    for line in pdb_lines:
        key = extract_atom_key_from_line(line, res_map)
        coord_string = line[30:54]  # keep x, y, z, exactly
        if key in coord_map:
            raise ValueError(f"Duplicate atom key found: {key}")
        coord_map[key] = coord_string
    return coord_map

def merge_coordinates(ref_pdb, pdb, output_file):
    """Merge coordinates based on atom name and residue mapping."""
    with open(ref_pdb) as f:
        ref_lines = f.readlines()
    with open(pdb) as f:
        pdb_lines = f.readlines()

    ref_atoms = [line for line in ref_lines if line.startswith(("ATOM", "HETATM"))]
    pdb_atoms = [line for line in pdb_lines if line.startswith(("ATOM", "HETATM"))]

    # Build residue maps and coordinate map
    ref_res_map = build_residue_renumbering_map(ref_atoms)
    pdb_res_map = build_residue_renumbering_map(pdb_atoms)
    coord_map = build_coord_line_map(pdb_atoms, pdb_res_map)

    updated_lines = []
    for line in ref_lines:
        if line.startswith(('ATOM', 'HETATM')):
            key = extract_atom_key_from_line(line, ref_res_map)
            if key not in coord_map:
                raise ValueError(f"❌ Atom not found in coordinate target PDB: {key}")
            coord_str = coord_map[key]
            new_line = line[:30] + coord_str + line[54:]
            updated_lines.append(new_line)
        else:
            updated_lines.append(line)

    with open(output_file, "w") as out:
        out.writelines(updated_lines)

    print(f"✅ Coordinates updated and saved to '{output_file}'.")

def main():
    parser = argparse.ArgumentParser(
        description="Replace coordinates in a reference PDB using atom mapping (atom name, residue name, chain, relative index)."
    )
    parser.add_argument("--ref", required=True, help="Reference PDB file with correct atom info.")
    parser.add_argument("input", help="New PDB file with coordinates to use.")
    parser.add_argument("-o", "--output", default="merged.pdb", help="Output PDB file (default: merged.pdb)")
    args = parser.parse_args()

    merge_coordinates(args.ref, args.input, args.output)

if __name__ == "__main__":
    main()

