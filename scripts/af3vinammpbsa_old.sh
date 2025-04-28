#!/bin/bash

START_TIME=$SECONDS  # Start timer

PDBID_LigandID=$1

scripts_dir="/home/xg69107/program/pfas_docking/scripts"
af3_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID"
vina_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID/input/best_pose"

# Native
mkdir -p native
cp $af3_dir/input.pdb native_model.pdb
bash $scripts_dir/get_ref_for_af3vinammpbsa.sh --input_pdb native_model.pdb --work_dir native

# AF3
mkdir -p af3
cp $af3_dir/input.json af3.json
cp $af3_dir/input/input_model.cif af3/model.cif
mkdir -p af3/best_pose;       cp $af3_dir/input/best_pose/model.cif        af3/best_pose/model.cif
mkdir -p af3/seed-1_sample-0; cp $af3_dir/input/seed-1_sample-0/model.cif  af3/seed-1_sample-0/model.cif
mkdir -p af3/seed-1_sample-1; cp $af3_dir/input/seed-1_sample-1/model.cif  af3/seed-1_sample-1/model.cif
mkdir -p af3/seed-1_sample-2; cp $af3_dir/input/seed-1_sample-2/model.cif  af3/seed-1_sample-2/model.cif
mkdir -p af3/seed-1_sample-3; cp $af3_dir/input/seed-1_sample-3/model.cif  af3/seed-1_sample-3/model.cif
mkdir -p af3/seed-1_sample-4; cp $af3_dir/input/seed-1_sample-4/model.cif  af3/seed-1_sample-4/model.cif
bash $scripts_dir/af3.sh --input_json af3.json --native_dir native --run_af3 false

# Vina
mkdir -p vina
cp $af3_dir/x_ligand.pdb             vina/vina_ligand.pdb
cp $af3_dir/x_receptor.pdb           vina/vina_receptor.pdb
cp $vina_dir/x_ligand.pdbqt          vina/vina_ligand.pdbqt
cp $vina_dir/x_ligand_docked.pdbqt   vina/vina_ligand_docked.pdbqt
bash $scripts_dir/vina.sh --input_pdb native/native_model.pdb --work_dir vina --native_dir native --run_vina false
for ligand_pdb in $(ls vina/vina_ligand_*_convert.pdb | sort); do
    num=$(basename "$ligand_pdb" | sed -E 's/.*_ligand_([0-9]+)_convert\.pdb/\1/')
    cat vina/vina_receptor.pdb "$ligand_pdb" | grep 'ATOM' > vina/vina_model_${num}_convert.pdb
done

# MMPBSA
bash $scripts_dir/run_all_mmpbsa.sh

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') ✅ Total time taken: $ELAPSED_TIME seconds"

