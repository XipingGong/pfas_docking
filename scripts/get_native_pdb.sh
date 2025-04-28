#!/bin/bash


START_TIME=$SECONDS

# This bash script "get_native_pdb.sh" is used to create a native protein-ligand structure.
# This script has been successfully tested in the Sapelo2@GACRC.
# To run it, you may need to change the INPUT section in this script.

# INPUT
# -----
scripts_dir='/home/xg69107/program/pfas_docking/scripts'
data_dir='/home/xg69107/program/pfas_docking/data'
python="/home/xg69107/program/anaconda/anaconda3/bin/python" # it requires to install MDTraj, etc.

# INFO
# ----
print_help() {
    echo ""
    echo "Usage: bash get_native_pdb.sh [--pdbid PDBID] [--ligandid LIGANDID] [--work_dir DIR]"
    echo ""
    echo "Optional arguments:"
    echo "  --pdbid         Protein PDB ID (default: 7FEU)"
    echo "  --ligandid      Ligand 3-letter code (default: 4I6)"
    echo "  --work_dir      Working directory name (default: current directory)"
    echo "  -h, --help      Show this help message and exit"
    echo ""
    exit 0
}

# Default values
default_work_dir=$(pwd)
work_dir=$default_work_dir
PDBID="7FEU"
LigandID="4I6"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --pdbid) PDBID="$2"; shift ;;
        --ligandid) LigandID="$2"; shift ;;
        --work_dir) work_dir="$2"; shift ;;
        -h|--help) print_help ;;
        *) echo "❌ Unknown parameter: $1"; print_help ;;
    esac
    shift
done

PDBID_LigandID="${PDBID}_${LigandID}"

# Convert work_dir to full path
[[ "$work_dir" != /* ]] && work_dir="$default_work_dir/$work_dir"

# Validate $work_dir
if [ ! -d "$work_dir" ]; then
    echo "❌ Error: work_dir does not exist: $work_dir"
    exit 1
fi

output_pdb="${PDBID}_${LigandID}.pdb"
output_pdb_base="${output_pdb##*/}"
output_pdb_basename="${output_pdb_base%.*}"



# OPTIONAL (you may modify this section, but do not have to, except you know what you are doing)
# --------

# Starting
echo "# Creating a native structure >>"
echo ""

# Go to the $work_dir directory
echo "# Go to the $work_dir directory"
echo "$ cd $work_dir"
        cd $work_dir
echo ""

if [ -d "$work_dir" ] && [ "$(ls -A "$work_dir")" ]; then
    echo "⚠️  Warning: Working directory already exists and is not empty:"
    echo "   → $work_dir"
    echo "   Existing data may be overwritten."
else
    echo "✅ All data will be saved into this working directory:"
    echo "   - $work_dir"
fi
echo ""

echo "# Download a PDB file from RCSB PDB"
echo "$ $python $scripts_dir/download_pdb.py $PDBID --ligand_id $LigandID"
        $python $scripts_dir/download_pdb.py $PDBID --ligand_id $LigandID # "download_pdb.py" can be found in the "scripts" folder
echo ""

echo "# Save the downloaded pdb as ${PDBID_LigandID}_ori.pdb"
echo "mv ${PDBID_LigandID}.pdb ${PDBID_LigandID}_ori.pdb"
      mv ${PDBID_LigandID}.pdb ${PDBID_LigandID}_ori.pdb

echo "# Check what this PDB file has"
echo "$ $python $scripts_dir/check_pdb.py ${PDBID_LigandID}_ori.pdb"
        $python $scripts_dir/check_pdb.py ${PDBID_LigandID}_ori.pdb # check what this PDB file has
echo ""

# Clean up this PDB file, so that it only has one protein and one ligand molecule
echo "# Clean up this PDB file: ${PDBID_LigandID}_ori.pdb"
echo "------------------------"
echo ""

echo "# Fill missing heavy atoms in a PDB file"
echo "$ $python $scripts_dir/clean_pdb.py ${PDBID_LigandID}_ori.pdb"
        $python $scripts_dir/clean_pdb.py ${PDBID_LigandID}_ori.pdb # Fill missing heavy atoms in a PDB file
echo ""

echo "# Extract the protein and ligand"
echo "$ $python $scripts_dir/extract_pdb.py cleaned_${PDBID_LigandID}_ori.pdb $LigandID"
        $python $scripts_dir/extract_pdb.py cleaned_${PDBID_LigandID}_ori.pdb $LigandID # extract the protein and ligand
echo ""

echo "# Save it as the pdb file(s)"
pdb_files=( $(ls cleaned_*_ori_*.pdb 2>/dev/null | sort) )

n=${#pdb_files[@]}
if (( n == 0 )); then
    echo "❌ No cleaned PDB files found."
    exit 1
fi

for i in "${!pdb_files[@]}"; do
    old="${pdb_files[$i]}"
    if (( n == 1 )); then
        new="${output_pdb_basename}.pdb"
    else
        new="${output_pdb_basename}_$((i + 1)).pdb"
    fi

    echo "Renaming $old -> $new"
    mv "$old" "$new"
    echo ""

    echo "# Check this pdb file (it must have two chains)"
    $python "$scripts_dir/check_pdb.py" "$new"
    echo ""
done

#>> native_model.pdb: <mdtraj.Trajectory with 1 frames, 1072 atoms, 134 residues, and unitcells>
#>>
#>> Protein Info:
#>>   - Chain 'A':
#>>     - Number of residues: 133
#>>     - Number of atoms: 1044
#>>   - native_model.pdb: Multiple protein chains present: False
#>>
#>> Ligand Info (including water):
#>>   - Ligand '4I6':
#>>     - Chain IDs: ['B']
#>>     - Number of residues: 1
#>>     - Number of atoms: 28

echo "# Remove the temporary files"
echo "$ rm -f \#* x_* cleaned_*"
        rm -f \#* x_* cleaned_*
echo ""

echo "# Go back to the current directory"
echo "$ cd $default_work_dir"
        cd $default_work_dir # go back to the current directory
echo ""

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for get_native_pdb.sh: $ELAPSED_TIME seconds"

