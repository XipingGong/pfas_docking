#!/bin/bash


START_TIME=$SECONDS
python="/home/xg69107/program/anaconda/anaconda3/bin/python"
scripts_dir='/home/xg69107/program/pfas_docking/scripts' # need to change 

mkdir -p mmpbsa

# Function to process each model
process_model() {
    local pdb_file="$1"
    local base_name="${pdb_file%.pdb}"
    local work_pdb="mmpbsa/$pdb_file"

    echo ">> Processing: $pdb_file"

    # Step 1: Generate H-added input pdb file
    bash $scripts_dir/get_ref_for_af3vinammpbsa.sh \
        --input_pdb "$work_pdb" \
        --work_dir mmpbsa \
        --native_dir native

    # Step 2: Run MMPBSA
    bash $scripts_dir/mmpbsa.sh \
        --input_pdb "mmpbsa/${base_name}H.pdb" \
        --work_dir mmpbsa \
        --native_dir native

    # Step 3: Check RMSD values
    if [[ $base_name == af3* ]]; then
        ref_file="mmpbsa/af3-best_pose-aligned_model_convertH_emin.pdb"
    elif [[ $base_name == vina* ]]; then
        ref_file="mmpbsa/vina-vina_model_1_convertH_emin.pdb"
    elif [[ $base_name == nat* ]]; then
        ref_file="native/native_modelH.pdb"
    else
        echo "⚠️ Warning: Unknown prefix for $base_name — skipping RMSD"
        return
    fi
    $python $scripts_dir/check_rmsd.py \
        --ref "$ref_file" \
        "mmpbsa/${base_name}H_emin.pdb" \
        > "mmpbsa/${base_name}H_emin_Top1RMSD.dat"

    # Step 4: Output MMPBSA & RMSD result
    echo "📄 Results for $base_name:"
    cat "mmpbsa/${base_name}H_emin_MMPBSA.dat"
    cat "mmpbsa/${base_name}H_emin_Top1RMSD.dat"
    echo "✅ Done: $pdb_file"
    echo "--------------------------------------------"
}

echo "🚀 Starting AF3 and Vina models..."
for pdb_file in {native/native_model.pdb,af3/*/aligned_model_convert.pdb,vina/vina_model_*_convert.pdb}; do
    if [[ -f "$pdb_file" ]]; then
        full_path=$(realpath "$pdb_file")
        safe_name=$(echo "$pdb_file" | sed 's|/|-|g')

        cp "$full_path" "mmpbsa/$safe_name"
        process_model "$safe_name"

    fi
done

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for run_all_mmpbsa.sh script: $ELAPSED_TIME seconds"

