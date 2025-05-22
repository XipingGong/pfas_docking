#!/bin/bash

START_TIME=$SECONDS

# This bash script "mmpbsa.sh" is used to run gmx_MM/PBSA for calculating protein-ligand binding free energy using a single trajectory (ST) protocol.
# This script has been successfully tested in the computing center @GACRC.
# To run it, you may need to change the INPUT section in this script.

# INPUT
# -----

# This is a "scripts" folder that can be downloaded from https://github.com/XipingGong/pfas_docking.git
scripts_dir='/home/xg69107/program/pfas_docking/scripts' # need to change 

# These should be installed
python="/home/xg69107/program/anaconda/anaconda3/bin/python"
gmx="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin.AVX2_256/gmx"
gmx_MMPBSA="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA" # need to run in the gmxMMPBSA environment

# INFO
# ----
print_help() {
    echo ""
    echo "Usage: bash mmpbsa.sh --input_pdb FILE [--work_dir DIR] [--native_dir DIR]"
    echo ""
    echo "Required:"
    echo "  --input_pdb     PDB file containing protein (chain A) and ligand (chain B), with hydrogens"
    echo ""
    echo "Optional:"
    echo "  --work_dir      Working directory (default: current directory)"
    echo "  --native_dir    A directory to include the files as the reference (default: native in the current directory)"
    echo "  -h, --help      Show this help message and exit"
    echo ""
    exit 0
}

# Default Values
default_work_dir=$(pwd)
work_dir=$default_work_dir

# This pdb file has one protein and one ligand.
# The file name must be in the working directory.
# It must have two chains, like chain A for protein and chain B for ligand.
input_pdb=""

# Parse Arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pdb) input_pdb="$2"; shift ;;
        --work_dir) work_dir="$2"; shift ;;
        --native_dir) native_dir="$2"; shift ;;
        -h|--help) print_help ;;
        *) echo "❌ Unknown parameter: $1"; print_help ;;
    esac
    shift
done
if [ -z "$input_pdb" ]; then
    echo "❌ Error: --input_pdb is required."
    print_help
fi

# Convert input_pdb to full path
[[ "$input_pdb" != /* ]] && input_pdb="$PWD/$input_pdb"

# Validate $input_pdb
if [ ! -f "$input_pdb" ]; then
    echo "❌ Error: $input_pdb must be a valid file."
    print_help
fi
input_pdb_base="${input_pdb##*/}"
input_pdb_basename="${input_pdb_base%.*}"

# Convert work_dir to full path
[[ "$work_dir" != /* ]] && work_pdb="$PWD/$work_dir"

# Validate $work_dir
if [ ! -d "$work_dir" ]; then
    echo "❌ Error: work_dir does not exist: $work_dir"
    exit 1
fi

# Convert native_dir to full path
[[ "$native_dir" != /* ]] && native_dir="$PWD/$native_dir"

# Validate $native_dir
if [ ! -d "$native_dir" ]; then
    echo "❌ Error: native_dir does not exist: $native_dir"
    exit 1
fi

# OPTIONAL (you could need to modify Vina Section, but do not modify others, except you know what you are doing)
# --------

# Starting
echo "# Running a vaccum energy minimization job >>"
echo ""

# Go to the $work_dir directory
echo "$ cd $work_dir"
        cd $work_dir
echo ""

# Vacuum energy minimization
# --------------------------

echo "$gmx editconf -f $input_pdb -o complex_init_box.pdb -bt cubic -d 5.0 "
      $gmx editconf -f $input_pdb -o complex_init_box.pdb -bt cubic -d 5.0 
echo ""

# steep-optimization
echo "$gmx grompp -f $native_dir/em_steep.mdp -c complex_init_box.pdb -r complex_init_box.pdb -p $native_dir/native_modelH.top -o em_steep.tpr -maxwarn 1"
      $gmx grompp -f $native_dir/em_steep.mdp -c complex_init_box.pdb -r complex_init_box.pdb -p $native_dir/native_modelH.top -o em_steep.tpr -maxwarn 1
echo ""

echo "$gmx mdrun -v -deffnm em_steep"
      $gmx mdrun -v -deffnm em_steep
echo ""

# cg-optimization
echo "$gmx grompp -f $native_dir/em_cg.mdp -c em_steep.gro -r complex_init_box.pdb -p $native_dir/native_modelH.top -o em_cg.tpr -maxwarn 1"
      $gmx grompp -f $native_dir/em_cg.mdp -c em_steep.gro -r complex_init_box.pdb -p $native_dir/native_modelH.top -o em_cg.tpr -maxwarn 1
echo ""

echo "$gmx mdrun -v -deffnm em_cg"
      $gmx mdrun -v -deffnm em_cg
echo ""

echo "echo -e "0\n0" | $gmx trjconv -f em_cg.gro -s em_cg.tpr -o em_cg_mol.pdb -pbc mol -center -ur compact -dump 0"
      echo -e "0\n0" | $gmx trjconv -f em_cg.gro -s em_cg.tpr -o em_cg_mol.pdb -pbc mol -center -ur compact -dump 0
echo ""

# - Aligned complex_emin.pdb
echo "$python $scripts_dir/align_pdb.py em_cg_mol.pdb --ref $native_dir/native_modelH.pdb --output aligned_modelH.pdb --oligand aligned_ligandH.pdb"
      $python $scripts_dir/align_pdb.py em_cg_mol.pdb --ref $native_dir/native_modelH.pdb --output aligned_modelH.pdb --oligand aligned_ligandH.pdb
echo ""

echo "$python $scripts_dir/update_pdb_coord.py aligned_modelH.pdb --ref $native_dir/native_modelH.pdb -o ${input_pdb_basename}_emin.pdb"
      $python $scripts_dir/update_pdb_coord.py aligned_modelH.pdb --ref $native_dir/native_modelH.pdb -o ${input_pdb_basename}_emin.pdb
echo ""

# MM/PBSA
echo "# Run the MM/PBSA (change it if needed)"
echo "echo -e "name 1 receptor\nname 13 ligand\n1 | 13\nq" | $gmx make_ndx -f em_cg.tpr -o index.ndx"
      echo -e "name 1 receptor\nname 13 ligand\n1 | 13\nq" | $gmx make_ndx -f em_cg.tpr -o index.ndx
echo ""
echo "$gmx_MMPBSA -O -i $native_dir/mmpbsa.in -cs em_cg.tpr -ct em_cg.trr -ci index.ndx -cg 1 13 -cp $native_dir/native_modelH.top -nogui -o ${input_pdb_basename}_emin_MMPBSA.dat -eo ${input_pdb_basename}_emin_MMPBSA.csv"
      $gmx_MMPBSA -O -i $native_dir/mmpbsa.in -cs em_cg.tpr -ct em_cg.trr -ci index.ndx -cg 1 13 -cp $native_dir/native_modelH.top -nogui -o ${input_pdb_basename}_emin_MMPBSA.dat -eo ${input_pdb_basename}_emin_MMPBSA.csv
echo ""
echo "cat ${input_pdb_basename}_emin_MMPBSA.dat"
      cat ${input_pdb_basename}_emin_MMPBSA.dat
echo ""

# Analysis
# --------

# Check RMSD values
echo "$python "$scripts_dir/check_rmsd.py" --ref $native_dir/native_modelH.pdb "${input_pdb_basename}_emin.pdb" | tee "${input_pdb_basename}_emin_RMSD.dat""
      $python "$scripts_dir/check_rmsd.py" --ref $native_dir/native_modelH.pdb "${input_pdb_basename}_emin.pdb" | tee "${input_pdb_basename}_emin_RMSD.dat"
echo ""

echo "# Remove the temporary files"
echo "$ rm -f \#* x_* "
        rm -f \#* x_* 
echo ""

echo "# Go back to the current directory"
echo "$ cd $default_work_dir"
        cd $default_work_dir
echo ""

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for mmpbsa.sh: $ELAPSED_TIME seconds"

