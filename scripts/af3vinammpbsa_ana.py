#!/usr/bin/env python3

import argparse
import re
import sys
from collections import defaultdict

def parse_blocks(input_file):
    blocks = defaultdict(list)
    current_block = None
    has_processing = False

    with open(input_file, 'r') as f:
        for line in f:
            if "Processing" in line:
                has_processing = True
                match = re.search(r'Processing (\S+)', line)
                if match:
                    current_block = match.group(1)
            if current_block:
                blocks[current_block].append(line.strip())

    # If no 'Processing' found, treat the whole file as one block
    if not has_processing:
        blocks["SingleExample"] = []
        with open(input_file, 'r') as f:
            for line in f:
                blocks["SingleExample"].append(line.strip())

    return blocks

def process_block(block_lines):
    top1rmsd_data = {}
    mmpbsa_data = {}
    rmsd_data = {}

    for line in block_lines:
        if re.search(r'_MMPBSA\.dat:ΔTOTAL', line):
            file = line.split(':')[0]
            energy = float(line.split(':')[1].split()[1])
            mmpbsa_data[file] = energy

        elif re.search(r'_Top1RMSD\.dat:.*Ligand RMSD \(Direct\)', line):
            file = line.split(':')[0]
            rmsd = float(re.search(r'\[([0-9\.]*)\]', line).group(1))
            top1rmsd_data[file] = rmsd

        elif re.search(r'_RMSD\.dat:.*Ligand RMSD \(Direct\)', line):
            file = line.split(':')[0]
            rmsd = float(re.search(r'\[([0-9\.]*)\]', line).group(1))
            rmsd_data[file] = rmsd

    return top1rmsd_data, mmpbsa_data, rmsd_data

def analyze_block(name, top1rmsd_data, mmpbsa_data, rmsd_data):
    qualified_files = [f for f, rmsd in top1rmsd_data.items() if rmsd < 0.2]

    if not qualified_files:
        print(f"❌ {name}: No pose with Top1RMSD < 0.2 nm.\n")
        return

    best_file = min(qualified_files, key=lambda f: mmpbsa_data.get(f.replace('_Top1RMSD', '_MMPBSA'), float('inf')))

    mmpbsa_file = best_file.replace('_Top1RMSD', '_MMPBSA')
    rmsd_file = best_file.replace('_Top1RMSD', '_RMSD')

    mmpbsa_value = mmpbsa_data.get(mmpbsa_file, None)
    top1rmsd_value = top1rmsd_data.get(best_file, None)
    rmsd_value = rmsd_data.get(rmsd_file, None)

    print(f"✅ Best pose (AF3-Vina-MMPBSA) - {best_file[:-13]}.pdb ; MMPBSA ΔTOTAL: {mmpbsa_value:.2f} kcal/mol; Top1RMSD: {top1rmsd_value:.3f} nm ; RMSD: {rmsd_value:.3f} nm\n")

def main():
    parser = argparse.ArgumentParser(description="Analyze docking results per molecule: find best pose with Top1RMSD < 0.2 nm and minimum ΔTOTAL energy.")
    parser.add_argument("input_file", help="Input log file containing one or multiple Processing blocks")
    args = parser.parse_args()

    blocks = parse_blocks(args.input_file)

    if not blocks:
        print("❌ No valid data found in the input file.")
        sys.exit(1)

    for name, lines in blocks.items():
        top1rmsd_data, mmpbsa_data, rmsd_data = process_block(lines)
        analyze_block(name, top1rmsd_data, mmpbsa_data, rmsd_data)

if __name__ == "__main__":
    main()

