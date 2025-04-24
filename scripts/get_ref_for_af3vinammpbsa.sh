#!/bin/bash

START_TIME=$SECONDS

# This bash script "get_ref_for_af3vinammpbsa.sh" is used to create the input files as the reference for running an AF3-Vina-MMPBSA job.
# This script has been successfully tested in the Sapelo2@GACRC.
# To run it, you may need to change the INPUT section in this script.

# INPUT
# -----

scripts_dir='/home/xg69107/program/pfas_docking/scripts'
data_dir='/home/xg69107/program/pfas_docking/data'
export GMXLIB="$data_dir/gmxff"
protein_ff="amber14sb_OL15"

# These should be installed
gmx="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin.AVX2_256/gmx"
acpype="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin/acpype" # run it in the gmxMMPBSA environment
obabel="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin/obabel"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"


# INFO
# ----
print_help() {
    echo ""
    echo "Usage: bash get_ref_for_af3vinammpbsa.sh --input_pdb FILE [--work_dir DIR] [--native_dir DIR]"
    echo ""
    echo "Required:"
    echo "  --input_pdb     PDB file containing protein (chain A) and ligand (chain B), without hydrogens"
    echo ""
    echo "Optional arguments:"
    echo "  --work_dir      Working directory (default: current directory)"
    echo "  --native_dir    Directory to save the reference structures (default: work_dir)"
    echo "                  Note: if the native_dir exists, then it will add hydrogens and a new input_pdbH.pdb will be created"
    echo "  -h, --help      Show this help message and exit"
    echo ""
    exit 0
}

# Default Values
default_work_dir=$(pwd)
work_dir=$default_work_dir
native_dir=$work_dir
has_native_dir="false"

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

# Check $input_pdb
if [ -z "$input_pdb" ]; then
    echo "❌ Error: --input_pdb is required."
    print_help
fi

# Convert $input_pdb to full path
[[ "$input_pdb" != /* ]] && input_pdb="$PWD/$input_pdb"

# Validate $input_pdb
if [ ! -f "$input_pdb" ]; then
    echo "❌ Error: $input_pdb must be a valid file."
    print_help
fi
input_pdb_base="${input_pdb##*/}" # e.g., vina.pdb
input_pdb_basename="${input_pdb_base%.*}" # e.g., vina

# Convert work_dir to full path
[[ "$work_dir" != /* ]] && work_dir="$default_work_dir/$work_dir"

# Validate $work_dir
if [ ! -d "$work_dir" ]; then
    echo "❌ Error: work_dir does not exist: $work_dir"
    exit 1
fi

# Convert native_dir to full path
[[ "$native_dir" != /* ]] && native_dir="$default_work_dir/$native_dir"

# Validate $native_dir
if [ ! -d "$native_dir" ]; then
    echo "❌ Error: native_dir does not exist: $native_dir"
    exit 1
fi

# OPTIONAL (you may modify this section, but do not have to, except you know what you are doing)
# --------

#
# Starting
echo "# Creating the reference structures >>"
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

native_proteinH_pdb="$native_dir/native_proteinH.pdb"
native_proteinH_mol2="$native_dir/native_proteinH.mol2"
native_ligandH_pdb="$native_dir/native_ligandH.pdb"
native_ligandH_mol2="$native_dir/native_ligandH.mol2"

# Check if both files exist in the native directory
if [[ -f "$native_ligandH_pdb" && -f "$native_ligandH_mol2" && -f "$native_proteinH_pdb" && -f "$native_proteinH_mol2" ]]; then
    has_native_dir=true
else
    has_native_dir=false
fi

# Case 1: create a hydrogen-added structure
if [ "$has_native_dir" == "true" ]; then
    echo "# Case 1: Create a hydrogen-added structure"
    # Protein
    echo "$ grep ' A ' $input_pdb > x_prot.pdb"
            grep ' A ' $input_pdb > x_prot.pdb # Native protein
    echo ""

    echo "$ $gmx pdb2gmx -f x_prot.pdb -o x_protH_gmx.pdb -p x_protH_gmx.top -ff $protein_ff -water tip3p"
            $gmx pdb2gmx -f x_prot.pdb -o x_protH_gmx.pdb -p x_protH_gmx.top -ff $protein_ff -water tip3p
    echo ""

    echo "$ $python $scripts_dir/addh.py x_protH_gmx.pdb --ref_pdb $native_proteinH_pdb --ref_mol2 $native_proteinH_mol2 -o x_protH.pdb"
            $python $scripts_dir/addh.py x_protH_gmx.pdb --ref_pdb $native_proteinH_pdb --ref_mol2 $native_proteinH_mol2 -o x_protH.pdb
    echo ""

    # Ligand
    echo "$ grep ' B ' $input_pdb > x_lig.pdb"
            grep ' B ' $input_pdb > x_lig.pdb # Native ligand
    echo ""

    echo "$ $python $scripts_dir/addh.py x_lig.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o x_ligH.pdb"
            $python $scripts_dir/addh.py x_lig.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o x_ligH.pdb
    echo ""

    # Protein-Ligand
    echo "$python $scripts_dir/merge_2pdbs.py x_protH.pdb x_ligH.pdb -o x_modelH.pdb"
          $python $scripts_dir/merge_2pdbs.py x_protH.pdb x_ligH.pdb -o x_modelH.pdb
    echo ""

    echo "$gmx editconf -f x_modelH.pdb -o ${input_pdb_basename}H.pdb -bt cubic -d 2.0 -noc"
          $gmx editconf -f x_modelH.pdb -o ${input_pdb_basename}H.pdb -bt cubic -d 2.0 -noc
    echo ""

fi

# Case 2: creating the native reference structures
if [ "$has_native_dir" == "false" ]; then
    echo "# Case 2: Creating the native reference structures"
    echo "cp $input_pdb native_model.pdb"
          cp $input_pdb native_model.pdb
    echo ""

    # Extract the protein and ligand
    echo "$ grep ' A ' native_model.pdb > native_protein.pdb"
            grep ' A ' native_model.pdb > native_protein.pdb # Native protein
    echo ""
    echo "$ grep ' B ' native_model.pdb > native_ligand.pdb"
            grep ' B ' native_model.pdb > native_ligand.pdb # Native ligand
    echo ""

    # Add hydrogen atoms and save the native protein and ligand structures as a reference
    # - Native Protein
    echo "$ $gmx pdb2gmx -f native_protein.pdb -o x_native_proteinH.pdb -p native_proteinH.top -ff $protein_ff -water tip3p"
            $gmx pdb2gmx -f native_protein.pdb -o x_native_proteinH.pdb -p native_proteinH.top -ff $protein_ff -water tip3p
    echo ""
    echo "$ $gmx editconf -f x_native_proteinH.pdb -o native_proteinH.pdb"
            $gmx editconf -f x_native_proteinH.pdb -o native_proteinH.pdb
    echo ""
    echo "$ $obabel -ipdb native_proteinH.pdb -omol2 -O native_proteinH.mol2"
            $obabel -ipdb native_proteinH.pdb -omol2 -O native_proteinH.mol2
    echo ""


    # - Native Ligand
    echo "$ $obabel native_ligand.pdb -O x_native_ligandH.pdb -p 7.4"
            $obabel native_ligand.pdb -O x_native_ligandH.pdb -p 7.4
    echo ""
    echo "$ $python $scripts_dir/merge_pdb.py native_ligand.pdb x_native_ligandH.pdb > native_ligandH.pdb"
            $python $scripts_dir/merge_pdb.py native_ligand.pdb x_native_ligandH.pdb > native_ligandH.pdb
    echo ""
    echo "$ $obabel -ipdb native_ligandH.pdb -omol2 -O native_ligandH.mol2"
            $obabel -ipdb native_ligandH.pdb -omol2 -O native_ligandH.mol2
    echo ""

    echo "# A quick ligand structural optimization using MMFF94"
    echo "$obabel -imol2 native_ligandH.mol2 -omol2 -O ligandH.mol2 --minimize --ff MMFF94"
          $obabel -imol2 native_ligandH.mol2 -omol2 -O ligandH.mol2 --minimize --ff MMFF94 # necessary for ACPYPE running
    echo ""

    echo "# Running ACPYPE to generate ligand topology..."
    ligandH_charge=$(awk '/@<TRIPOS>ATOM/{f=1;next} /@<TRIPOS>/{f=0} f && NF>=9 {s+=$9} END {printf "%d", s + 0.5*(s>0) - 0.5*(s<0)}' ligandH.mol2)
    echo "✅ Ligand net charge: $ligandH_charge"
    echo "$acpype -i ligandH.mol2 -c bcc -n $ligandH_charge"
          $acpype -i ligandH.mol2 -c bcc -n $ligandH_charge # AM1-BCC charge
    echo ""
   #echo "$acpype -i ligandH.mol2 -c gas -n $ligandH_charge"
   #      $acpype -i ligandH.mol2 -c gas -n $ligandH_charge # Gasteiger charge
   #echo ""
    cp ligandH.acpype/posre_ligandH.itp posre_ligandH.itp
    cp ligandH.acpype/ligandH_GMX.top ligandH_GMX.top
    cp ligandH.acpype/ligandH_GMX.itp ligandH_GMX.itp

    # - Native Protein-Ligand complex
    echo "$ $python $scripts_dir/merge_2pdbs.py native_proteinH.pdb native_ligandH.pdb -o x_native_modelH.pdb"
            $python $scripts_dir/merge_2pdbs.py native_proteinH.pdb native_ligandH.pdb -o x_native_modelH.pdb
    echo ""
    echo "$ $gmx editconf -f x_native_modelH.pdb -o native_modelH.pdb -bt cubic -d 2.0 -noc"
            $gmx editconf -f x_native_modelH.pdb -o native_modelH.pdb -bt cubic -d 2.0 -noc
    echo ""
    echo "$ $python $scripts_dir/merge_2tops.py native_proteinH.top ligandH_GMX.top > native_modelH.top"
            $python $scripts_dir/merge_2tops.py native_proteinH.top ligandH_GMX.top > native_modelH.top
    echo ""

    # Copy the additional required input files
    echo "# Copy the additional required input files"
    echo "cp $scripts_dir/em.mdp $work_dir/em.mdp"
          cp $scripts_dir/em.mdp $work_dir/em.mdp
    echo ""
    echo "cp $scripts_dir/mmpbsa.in $work_dir/mmpbsa.in"
          cp $scripts_dir/mmpbsa.in $work_dir/mmpbsa.in
    echo ""

fi

# Remove the temporary files
echo "$ rm -f \#* x_*"
        rm -f \#* x_*
echo ""

# Go back to the current directory
echo "# Go back to the current directory"
echo "$ cd $default_work_dir"
        cd $default_work_dir
echo ""

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for get_ref_for_af3vinammpbsa.sh: $ELAPSED_TIME seconds"
echo ""

