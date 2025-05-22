#!/bin/bash

START_TIME=$SECONDS  # Start timer

PDBID_LigandID=$1

scripts_dir="/home/xg69107/program/pfas_docking/scripts"
af3_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID"
vina_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID/input/best_pose"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"

# AF3
#mkdir -p af3
#cp $af3_dir/input.json af3.json
#cp $af3_dir/input/input_model.cif af3/model.cif
#mkdir -p af3/best_pose;       cp $af3_dir/input/best_pose/model.cif        af3/best_pose/model.cif
#mkdir -p af3/seed-1_sample-0; cp $af3_dir/input/seed-1_sample-0/model.cif  af3/seed-1_sample-0/model.cif
#mkdir -p af3/seed-1_sample-1; cp $af3_dir/input/seed-1_sample-1/model.cif  af3/seed-1_sample-1/model.cif
#mkdir -p af3/seed-1_sample-2; cp $af3_dir/input/seed-1_sample-2/model.cif  af3/seed-1_sample-2/model.cif
#mkdir -p af3/seed-1_sample-3; cp $af3_dir/input/seed-1_sample-3/model.cif  af3/seed-1_sample-3/model.cif
#mkdir -p af3/seed-1_sample-4; cp $af3_dir/input/seed-1_sample-4/model.cif  af3/seed-1_sample-4/model.cif
#bash $scripts_dir/af3.sh --input_json af3.json --native_dir native --run_af3 false

# Native
mkdir -p native; rm -rf native/*
cp $af3_dir/input.pdb native/input.pdb
cp "af3/seed-1_sample-4/aligned_model_convert.pdb" native/native_model.pdb
bash $scripts_dir/get_ref_for_af3vinammpbsa.sh --input_pdb native/native_model.pdb --work_dir native

# Vina
mkdir -p vina; rm -rf vina/*
bash $scripts_dir/vina.sh --input_pdb native/native_model.pdb --work_dir vina --native_dir native --run_vina true

label="af3_s1_s4"
# MMPBSA
mkdir -p mmpbsa; #rm -rf mmpbsa/*
for ligand_path in $(ls vina/vina_ligand_*_convert.pdb | sort); do
    ligand_pdb=$(echo "$ligand_path" | sed 's|/|-|g')
    ligand_pdb_basename="${ligand_pdb%.pdb}"
    input_pdb="${label}_vina_model_${ligand_pdb_basename}.pdb"
    input_pdb_basename="${input_pdb%.pdb}"
    cat "vina/vina_receptor.pdb" "$ligand_path" | grep 'ATOM' > "mmpbsa/$input_pdb"
    bash "$scripts_dir/get_ref_for_af3vinammpbsa.sh" --input_pdb "mmpbsa/$input_pdb" --work_dir mmpbsa --native_dir native
    bash "$scripts_dir/mmpbsa.sh" --input_pdb "mmpbsa/${input_pdb_basename}H.pdb" --work_dir mmpbsa --native_dir native
    $python "$scripts_dir/check_rmsd.py" --ref native/input.pdb "mmpbsa/${input_pdb_basename}H_emin.pdb" > "mmpbsa/${input_pdb_basename}H_emin_NATRMSD.dat"
done

for ligand_path in $(ls af3/*/aligned_ligand_convert.pdb | sort); do
    sed 's/^\(.\{21\}\)A/\1B/' $ligand_path > "mmpbsa/x_ligand.pdb"
    ligand_pdb=$(echo "$ligand_path" | sed 's|/|-|g')
    ligand_pdb_basename="${ligand_pdb%.pdb}"
    input_pdb="${label}_af3_model_${ligand_pdb_basename}.pdb"
    input_pdb_basename="${input_pdb%.pdb}"
    cat "vina/vina_receptor.pdb" "mmpbsa/x_ligand.pdb" | grep 'ATOM' > "mmpbsa/$input_pdb"
    bash "$scripts_dir/get_ref_for_af3vinammpbsa.sh" --input_pdb "mmpbsa/$input_pdb" --work_dir mmpbsa --native_dir native
    bash "$scripts_dir/mmpbsa.sh" --input_pdb "mmpbsa/${input_pdb_basename}H.pdb" --work_dir mmpbsa --native_dir native
    $python "$scripts_dir/check_rmsd.py" --ref native/input.pdb "mmpbsa/${input_pdb_basename}H_emin.pdb" > "mmpbsa/${input_pdb_basename}H_emin_NATRMSD.dat"
done

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') ✅ Total time taken: $ELAPSED_TIME seconds"

