#!/bin/bash

START_TIME=$SECONDS  # Start timer

PDBID_LigandID=$1

scripts_dir="/home/xg69107/program/pfas_docking/scripts"
af3_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID"
vina_dir="/home/xg69107/work/pfas_pdbs/dock_dir/$PDBID_LigandID/input/best_pose"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"

# Native
mkdir -p native; rm -rf native/*
cp $af3_dir/input.pdb native/input.pdb
cp "$af3_dir/input.pdb" native/native_model.pdb #cp "af3/best_pose/aligned_model_convert.pdb" native/native_model.pdb
bash $scripts_dir/get_ref_for_af3vinammpbsa.sh --input_pdb native/native_model.pdb --work_dir native

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
mkdir -p vina; rm -rf vina/*
awk 'substr($0, 22, 1) == "A"' native/native_modelH.pdb > vina/x_protH.pdb # Native protein
awk 'substr($0, 22, 1) == "B"' af3/best_pose/aligned_model_convert.pdb > vina/x_lig.pdb # AF3-predicted best ligand
cat vina/x_protH.pdb vina/x_lig.pdb > vina/nat_prot_af3_best_lig.pdb
bash $scripts_dir/vina.sh --input_pdb vina/nat_prot_af3_best_lig.pdb --work_dir vina --native_dir native --run_vina true

# MMPBSA
mkdir -p mmpbsa; #rm -rf mmpbsa/*
for ligand_path in $(ls vina/vina_ligand_*_convert.pdb | sort); do
    ligand_pdb=$(echo "$ligand_path" | sed 's|/|-|g')
    ligand_pdb_basename="${ligand_pdb%.pdb}"
    input_pdb="vina_model_${ligand_pdb_basename}.pdb"
    input_pdb_basename="${input_pdb%.pdb}"
    cat "vina/vina_receptor.pdb" "$ligand_path" | grep 'ATOM' > "mmpbsa/$input_pdb"
    bash "$scripts_dir/get_ref_for_af3vinammpbsa.sh" --input_pdb "mmpbsa/$input_pdb" --work_dir mmpbsa --native_dir native
    bash "$scripts_dir/mmpbsa.sh" --input_pdb "mmpbsa/${input_pdb_basename}H.pdb" --work_dir mmpbsa --native_dir native
    $python "$scripts_dir/check_rmsd.py" --ref native/input.pdb "mmpbsa/${input_pdb_basename}H_emin.pdb" > "mmpbsa/${input_pdb_basename}H_emin_NATRMSD.dat"
    echo ">> Binding free energy estimation"
    file="mmpbsa/${input_pdb_basename}H_emin_MMPBSA.dat"
    echo "$ grep -E 'TOTAL|ΔTOTAL' "$file""
            grep -E 'TOTAL|ΔTOTAL' "$file"
    sum=$(grep -E 'TOTAL|ΔTOTAL' "$file" | tail -2 | awk '{s += $2} END {printf "%.2f", s}')
    echo "$file:MMPBSA_SUM(kcal/mol)= $sum" >> $file
    echo "MMPBSA_Sum $sum" >> $file
    tail -2 $file
    echo ""
done

for ligand_path in $(ls af3/*/aligned_ligand_convert.pdb | sort); do
    sed 's/^\(.\{21\}\)A/\1B/' $ligand_path > "mmpbsa/x_ligand.pdb"
    ligand_pdb=$(echo "$ligand_path" | sed 's|/|-|g')
    ligand_pdb_basename="${ligand_pdb%.pdb}"
    input_pdb="af3_model_${ligand_pdb_basename}.pdb"
    input_pdb_basename="${input_pdb%.pdb}"
    cat "vina/vina_receptor.pdb" "mmpbsa/x_ligand.pdb" | grep 'ATOM' > "mmpbsa/$input_pdb"
    bash "$scripts_dir/get_ref_for_af3vinammpbsa.sh" --input_pdb "mmpbsa/$input_pdb" --work_dir mmpbsa --native_dir native
    bash "$scripts_dir/mmpbsa.sh" --input_pdb "mmpbsa/${input_pdb_basename}H.pdb" --work_dir mmpbsa --native_dir native
    $python "$scripts_dir/check_rmsd.py" --ref native/input.pdb "mmpbsa/${input_pdb_basename}H_emin.pdb" > "mmpbsa/${input_pdb_basename}H_emin_NATRMSD.dat"
    echo ">> Binding free energy estimation"
    file="mmpbsa/${input_pdb_basename}H_emin_MMPBSA.dat"
    echo "$ grep -E 'TOTAL|ΔTOTAL' "$file""
            grep -E 'TOTAL|ΔTOTAL' "$file"
    sum=$(grep -E 'TOTAL|ΔTOTAL' "$file" | tail -2 | awk '{s += $2} END {printf "%.2f", s}')
    echo "$file:MMPBSA_SUM(kcal/mol)= $sum" >> $file
    echo "MMPBSA_Sum $sum" >> $file
    tail -2 $file
    echo ""
done

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') ✅ Total time taken: $ELAPSED_TIME seconds"

