#!/usr/bin/env python3
import glob
import os
import sys
import numpy as np
import re

# Check usage
if len(sys.argv) != 2:
    print("Usage: python tmp.py \"[pattern]\"")
    sys.exit(1)

pattern = sys.argv[1]
pdb_files = glob.glob(pattern)
if not pdb_files:
    print(f"No files found matching pattern: {pattern}")
    sys.exit(1)

# Define keywords (label, suffix, keyword in file)
terms = [
    # RMSD terms
    ("Protein_Backbone_RMSD_(Direct)", "_RMSD.dat", "Protein Backbone RMSD (Direct)"),
    ("Protein_Backbone_RMSD_(MDTraj)", "_RMSD.dat", "Protein Backbone RMSD (MDTraj)"),
    ("ΔLigand_RMSD_(Direct)", "_RMSD.dat", "Ligand RMSD (Direct)"),
    ("ΔLigand_RMSD_(MDTraj)", "_RMSD.dat", "Ligand RMSD (MDTraj)"),
    ("Protein_Backbone_Pocket_RMSD_(Direct)", "_RMSD.dat", "Protein Backbone Pocket RMSD (Direct)"),
    ("Protein_Backbone_Pocket_RMSD_(MDTraj)", "_RMSD.dat", "Protein Backbone Pocket RMSD (MDTraj)"),

    # NATRMSD terms
    ("Native_Protein_Backbone_RMSD_(Direct)", "_NATRMSD.dat", "Protein Backbone RMSD (Direct)"),
    ("Native_Protein_Backbone_RMSD_(MDTraj)", "_NATRMSD.dat", "Protein Backbone RMSD (MDTraj)"),
    ("Native_ΔLigand_RMSD_(Direct)", "_NATRMSD.dat", "Ligand RMSD (Direct)"),
    ("Native_ΔLigand_RMSD_(MDTraj)", "_NATRMSD.dat", "Ligand RMSD (MDTraj)"),
    ("Native_Protein_Backbone_Pocket_RMSD_(Direct)", "_NATRMSD.dat", "Protein Backbone Pocket RMSD (Direct)"),
    ("Native_Protein_Backbone_Pocket_RMSD_(MDTraj)", "_NATRMSD.dat", "Protein Backbone Pocket RMSD (MDTraj)"),

    # Energy terms
    ("ΔVDWAALS",   "_MMPBSA.dat", "ΔVDWAALS"),
    ("ΔEEL",       "_MMPBSA.dat", "ΔEEL"),
    ("ΔEGB",       "_MMPBSA.dat", "ΔEGB"),
    ("ΔESURF",     "_MMPBSA.dat", "ΔESURF"),
    ("ΔGGAS",      "_MMPBSA.dat", "ΔGGAS"),
    ("ΔGSOLV",     "_MMPBSA.dat", "ΔGSOLV"),
    ("ΔTOTAL",     "_MMPBSA.dat", "ΔTOTAL"),
    ("MMPBSA_Sum", "_MMPBSA.dat", "MMPBSA_Sum"),

]

# Add calculated headers
header = ["pdb_file"] + [label for label, _, _ in terms if label != "ΔTOTAL"] + ["ΔTOT_VDW", "ΔTOT_ELE", "ΔTOTAL", "ST_Flag"]

# Print table
print("#" + " ".join(header))

# Parse files
for pdb_file in sorted(pdb_files):
    file_prefix = os.path.splitext(pdb_file)[0]
    row = [file_prefix+'.pdb']
    values = {}

    for label, suffix, keyword in terms:
        target_file = file_prefix + suffix
        value = "NA"
        if os.path.isfile(target_file):
            with open(target_file, "r", encoding="utf-8") as f:
                for line in f:
                    if keyword in line:
                        try:
                            if "RMSD" in keyword:
                                match = re.search(r'Min\s*=\s*([\d\.Ee+-]+)', line)
                                if match:
                                    value = match.group(1)
                            else:
                                parts = line.split()
                               #print(f"debug: {line}")
                                for i, part in enumerate(parts):
                                    if part in keyword and i + 1 <= len(parts):
                                        value = parts[i + 1]
                                        break
                        except:
                            value = "NA"
                        break
        values[label] = value

    # Add all values except ΔTOTAL
    for label, _, _ in terms:
        if label != "ΔTOTAL":
            row.append(values.get(label, "NA"))

    # ΔTOT_VDW = ΔVDWAALS + ΔESURF
    try:
        tot_vdw = float(values["ΔVDWAALS"]) + float(values["ΔESURF"])
        row.append(f"{tot_vdw:.2f}")
    except:
        row.append("NA")

    # ΔTOT_ELE = ΔEEL + ΔEGB
    try:
        tot_ele = float(values["ΔEEL"]) + float(values["ΔEGB"])
        row.append(f"{tot_ele:.2f}")
    except:
        row.append("NA")

    # Add ΔTOTAL
    row.append(values.get("ΔTOTAL", "NA"))

    # Compute ST_Flag
    cutoff = 0.2
    ligand_rmsd = values.get("ΔLigand_RMSD_(MDTraj)", "NA")
    backbone_rmsd = values.get("Protein_Backbone_RMSD_(MDTraj)", "NA")

    try:
        ligand_rmsd_val = float(ligand_rmsd)
    except:
        ligand_rmsd_val = None

    try:
        backbone_rmsd_val = float(backbone_rmsd)
    except:
        backbone_rmsd_val = None

    if ligand_rmsd_val is not None and backbone_rmsd_val is not None:
        if ligand_rmsd_val < cutoff and backbone_rmsd_val < cutoff:
            flag = "Good"
        elif ligand_rmsd_val < cutoff and backbone_rmsd_val >= cutoff:
            flag = "Bad_ProteinBackbone"
        elif ligand_rmsd_val >= cutoff and backbone_rmsd_val < cutoff:
            flag = "Bad_Ligand"
        else:
            flag = "Bad_ProteinBackbone_Ligand"
    else:
        flag = "NA"

    row.append(flag)

    print(" ".join(row))

