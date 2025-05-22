#!/bin/bash

START_TIME=$SECONDS

python="/home/xg69107/program/anaconda/anaconda3/bin/python"
scripts_dir="/home/xg69107/program/pfas_docking/scripts"

show_help() {
    echo "Usage: bash run_all_mmpbsa.sh \"[pattern]\""
    echo ""
    echo "This script runs MMPBSA analysis on a set of PDB files that match the given glob pattern."
    echo ""
    echo "Arguments:"
    echo "  [pattern]         A glob pattern for input PDB files (in quotes)."
    echo "                    Example: \"vina/vina_model_*_convert.pdb\""
    echo ""
    echo "Options:"
    echo "  -h, --help        Show this help message and exit."
    echo ""
    echo "Output:"
    echo "  Processed files will be saved under the ./mmpbsa directory."
    echo ""
    exit 0
}

# Check for help
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    show_help
fi

# Check usage
if [ "$#" -ne 1 ]; then
    echo "❌ Error: Please provide exactly one glob pattern as argument."
    echo "Use --help for more information."
    exit 1
fi

input_pattern="$1"
pdb_files=( $input_pattern )  # Expand the pattern into an array

if [ ${#pdb_files[@]} -eq 0 ]; then
    echo "❌ No files found matching pattern: $input_pattern"
    exit 1
fi

mkdir -p mmpbsa

process_model() {
    local pdb_file="$1"
    local base_name="${pdb_file%.pdb}"
    local work_pdb="mmpbsa/$pdb_file"

    echo ">> Processing: $pdb_file"

    bash "$scripts_dir/get_ref_for_af3vinammpbsa.sh" \
        --input_pdb "$work_pdb" \
        --work_dir mmpbsa \
        --native_dir native

    bash "$scripts_dir/mmpbsa.sh" \
        --input_pdb "mmpbsa/${base_name}H.pdb" \
        --work_dir mmpbsa \
        --native_dir native

    echo "✅ Done: $pdb_file"
}

echo "🚀 Starting MMPBSA for files matching: $input_pattern"

for pdb_path in "${pdb_files[@]}"; do
    if [[ -f "$pdb_path" ]]; then
        full_path=$(realpath "$pdb_path")
        safe_name=$(echo "$pdb_path" | sed 's|/|-|g')
        cp "$full_path" "mmpbsa/$safe_name"
        process_model "$safe_name"
    else
        echo "⚠️ Skipping (not a file): $pdb_path"
    fi
done

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') ⏱️ Total time: $ELAPSED_TIME seconds"

