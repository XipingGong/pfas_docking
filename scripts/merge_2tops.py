#!/usr/bin/env python3
import argparse
import os
import re

def extract_ligandH_info(ligandH_top_file):
    """Extract include and molecule directives from a ligandH topology file, using full paths."""
    ligandH_dir = os.path.dirname(os.path.abspath(ligandH_top_file))
    include_itp = os.path.join(ligandH_dir, 'ligandH_GMX.itp')
    posre_itp  = os.path.join(ligandH_dir, 'posre_ligandH.itp')

    ligandH_include = f'#include "{include_itp}"\n'
    ligandH_posres  = (
        '; Ligand position restraints\n'
        '#ifdef POSRES_LIG\n'
        f'#include "{posre_itp}"\n'
        '#endif\n'
    )

    ligandH_molecule = ''
    with open(ligandH_top_file, 'r') as file:
        for line in file:
            if re.match(r"\s*ligandH\s+\d+", line):
                ligandH_molecule = line
                break

    if not ligandH_molecule:
        raise ValueError("Could not find 'ligandH' molecule directive in ligandH topology.")

    return ligandH_include, ligandH_posres, ligandH_molecule

def modify_top_file(input_file, ligandH_top_file):
    # Determine directories
    top_dir = os.path.dirname(os.path.abspath(input_file))
    ligandH_include, ligandH_posres, ligandH_molecule = extract_ligandH_info(ligandH_top_file)

    # Path for generic posre.itp in topol_protein.top
    posre_generic = os.path.join(top_dir, 'posre.itp')

    with open(input_file, 'r') as file:
        lines = file.readlines()

    found_ligand_include = False
    found_molecules = False
    in_molecules_section = False
    molecules_section = []
    modified_lines = []

    for line in lines:
        # Replace generic posre.itp include with full path
        if line.strip().startswith('#include') and 'posre.itp' in line and 'posre_ligandH.itp' not in line:
            line = f'#include "{posre_generic}"\n'

        # Capture [ molecules ] section
        if re.match(r'\[\s*molecules\s*\]', line):
            in_molecules_section = True
            found_molecules = True
            molecules_section.append(line)
            continue

        if in_molecules_section:
            molecules_section.append(line)
            continue

        # Append the original include line
        if line.strip().startswith('#include'):
            modified_lines.append(line)
            # Immediately follow with ligandH include
            if not found_ligand_include:
                modified_lines.append(ligandH_include)
                found_ligand_include = True
            # After generic posre include, insert ligandH posres block
            if 'posre.itp' in line and 'posre_ligandH.itp' not in line:
                modified_lines.append(ligandH_posres)
            continue

        # All other lines
        modified_lines.append(line)

    # Append ligandH molecule directive at end of molecules section
    if found_molecules:
        molecules_section.append(ligandH_molecule)
        modified_lines.extend(molecules_section)

    return ''.join(modified_lines)

def main():
    parser = argparse.ArgumentParser(
        description='Modify a GROMACS .top file to include full-path posre and ligandH directives.'
    )
    parser.add_argument('input_file',
                        help='Path to the input topology file (e.g., topol_protein.top)')
    parser.add_argument('ligandH_top_file',
                        help='Path to the ligandH topology file (e.g., ligandH.acpype/ligandH_GMX.top)')
    args = parser.parse_args()

    modified_content = modify_top_file(args.input_file, args.ligandH_top_file)
    print(modified_content, end='')

if __name__ == '__main__':
    main()

