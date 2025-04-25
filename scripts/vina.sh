#!/bin/bash

START_TIME=$SECONDS

# This bash script "vina.sh" is used to run AutoDock Vina for protein-ligand docking.
# This script has been successfully tested in the computing center @GACRC.
# To run it, you may need to change the INPUT section in this script.

# INPUT
# -----

# This is a "scripts" folder that can be downloaded from https://github.com/XipingGong/pfas_docking.git
scripts_dir='/home/xg69107/program/pfas_docking/scripts' # need to change 
obabel="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin/obabel"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"
vina="/home/xg69107/program/anaconda/anaconda3/bin/vina"
mk_prepare_ligand="/home/xg69107/program/anaconda/anaconda3/bin/mk_prepare_ligand.py"
mk_prepare_receptor="/home/xg69107/program/anaconda/anaconda3/bin/mk_prepare_receptor.py"


# INFO
# ----
print_help() {
    echo ""
    echo "Usage: bash vina.sh --input_pdb FILE [--work_dir DIR] [--native_dir FOLDER]"
    echo "                    [--receptor FILE] [--ligand FILE] [--pocket_params FILE]"
    echo ""
    echo "Required:"
    echo "  --input_pdb        PDB file containing protein (chain A) and ligand (chain B), without hydrogens"
    echo ""                    
    echo "Optional:"           
    echo "  --work_dir         Working directory (default: current directory)"
    echo "  --native_dir       Folder that contains the reference structures, such as native_model.pdb, native_ligand.pdb, native_ligandH.pdb, and native_ligandH.mol2,"
    echo "                     (default: if not provided, then the input_pdb file will be as a reference for the RMSD calculations without alignments)"
    echo "  --receptor         Pre-generated receptor file (default: vina_receptor.pdb created from input_pdb)"
    echo "  --ligand           Pre-generated ligand file (default: vina_ligand.pdb created from input_pdb)"
    echo "  --pocket_params    Pre-generated vina config file (default: vina_pocket_params.txt created from input_pdb)"
    echo "  -h, --help         Show this help message and exit"
    echo ""
    exit 0
}

# Default Values
default_work_dir=$(pwd)
work_dir=$default_work_dir

# This pdb file has one protein and one ligand, and their residue name must be in the RCSB PDB.
# The file name must be in the working directory.
# It must have two chains, like chain A for protein and chain B for ligand.
input_pdb=""
receptor_pdb=""
ligand_pdb=""
pocket_params=""
native_dir=""

# Parse Arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pdb) input_pdb="$2"; shift ;;
        --work_dir) work_dir="$2"; shift ;;
        --receptor) receptor_pdb="$2"; shift ;;
        --ligand) ligand_pdb="$2"; shift ;;
        --pocket_params) pocket_params="$2"; shift ;;
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

# Convert input to full path
[[ "$input_pdb" != /* ]] && input_pdb="$default_work_dir/$input_pdb"
[[ "$work_dir" != /* ]] && work_dir="$default_work_dir/$work_dir"
[[ "$native_dir" != /* && -n "$native_dir" ]] && native_dir="$default_work_dir/$native_dir"
[[ "$receptor_pdb" != /* ]] && receptor_pdb="$default_work_dir/$receptor_pdb"
[[ "$ligand_pdb" != /* ]] && ligand_pdb="$default_work_dir/$ligand_pdb"
[[ "$pocket_params" != /* ]] && pocket_params="$default_work_dir/$pocket_params"

# Validate $input_pdb
if [ ! -f "$input_pdb" ]; then
    echo "❌ Error: $input_pdb must be a valid file."
    print_help
fi
input_pdb_base="${input_pdb##*/}"
input_pdb_basename="${input_pdb_base%.*}"

# Validate $work_dir
if [ ! -d "$work_dir" ]; then
    echo "❌ Error: work_dir does not exist: $work_dir"
    exit 1
fi


# OPTIONAL (you could need to modify Vina Section, but do not modify others, except you know what you are doing)
# --------

# Starting
echo "# Running an AutoDock Vina job >>"
echo ""

# Create vina_dir for saving the data
vina_dir="$work_dir"

echo "# Go to the $vina_dir directory"
echo "$ cd $vina_dir"
        cd $vina_dir
echo ""

# Validate $vina_dir
if [ -d "$vina_dir" ] && [ "$(ls -A "$vina_dir")" ]; then
    echo "⚠️  Warning: Vina folder already exists and is not empty:"
    echo "   → $vina_dir"
    echo "   Existing data may be overwritten."
else
    echo "✅ All data will be saved into this folder:"
    echo "   → $vina_dir"
fi

echo "# Create the input files for Vina"
# Generate input files if not provided
if [[ -z "$receptor_pdb" || -z "$ligand_pdb" || -z "$pocket_params" ]]; then
    echo "$ $python $scripts_dir/get_inputs_for_vina.py $input_pdb --protein_output vina_receptor.pdb --ligand_output vina_ligand.pdb --pocket_params_output vina_pocket_params.txt"
            $python $scripts_dir/get_inputs_for_vina.py $input_pdb --protein_output vina_receptor.pdb --ligand_output vina_ligand.pdb --pocket_params_output vina_pocket_params.txt
    receptor_pdb="vina_receptor.pdb"
    ligand_pdb="vina_ligand.pdb"
    pocket_params="vina_pocket_params.txt"
else
    echo "> Using provided receptor, ligand, and pocket_params files:"
    echo "  - receptor: $receptor_pdb"
    echo "  - ligand: $ligand_pdb"
    echo "  - pocket_params: $pocket_params"
fi
echo ""

echo "$ $mk_prepare_receptor -i $receptor_pdb -o vina_receptor -p -a"
        $mk_prepare_receptor -i $receptor_pdb -o vina_receptor -p -a
echo ""
echo "$ $obabel $ligand_pdb -O vina_ligand.sdf -p 7.4"
        $obabel $ligand_pdb -O vina_ligand.sdf -p 7.4
echo ""
echo "$ $mk_prepare_ligand -i vina_ligand.sdf -o vina_ligand.pdbqt"
        $mk_prepare_ligand -i vina_ligand.sdf -o vina_ligand.pdbqt
echo ""

echo "# Run the Vina (change it if needed)"
echo "$ $vina --receptor vina_receptor.pdbqt --ligand vina_ligand.pdbqt --config vina_pocket_params.txt --out vina_ligand_docked.pdbqt --exhaustiveness 32 --num_modes 5"
        $vina --receptor vina_receptor.pdbqt \
              --ligand vina_ligand.pdbqt \
              --config $pocket_params \
              --out vina_ligand_docked.pdbqt \
              --exhaustiveness 32 \
              --num_modes 5
echo ""

# Analysis
# --------

echo "# Convert docked pdbqt to vina_ligands.pdb & individual vina_ligand_*.pdb"
echo "$ $obabel vina_ligand_docked.pdbqt -O vina_ligands.pdb --separate -d"
        $obabel vina_ligand_docked.pdbqt -O vina_ligands.pdb --separate -d
echo ""
echo "$ $obabel vina_ligands.pdb -O vina_ligand_.pdb -m -d"
        $obabel vina_ligands.pdb -O vina_ligand_.pdb -m -d
echo ""

echo "# Check their RMSD values: vina_ligands.pdb "
echo "$ $obabel vina_ligand.pdbqt -O vina_ligand_ref.pdb -d"
        $obabel vina_ligand.pdbqt -O vina_ligand_ref.pdb -d
echo ""
echo "$ $python $scripts_dir/check_rmsd.py --ref vina_ligand_ref.pdb vina_ligands.pdb"
        $python $scripts_dir/check_rmsd.py --ref vina_ligand_ref.pdb vina_ligands.pdb
echo ""

# Align the structures with the reference
# Why this step: the predicted structure files could miss the atom info included in the reference files.
native_model_pdb="$native_dir/native_model.pdb"
native_ligand_pdb="$native_dir/native_ligand.pdb"
native_ligandH_pdb="$native_dir/native_ligandH.pdb"
native_ligandH_mol2="$native_dir/native_ligandH.mol2"

if [[ -f "$native_model_pdb" && -f "$native_ligand_pdb" && -f "$native_ligandH_pdb" && -f "$native_ligandH_mol2" ]]; then

    echo "# ✅ All native reference files found. Proceeding with ligand alignment and adding hydrogen atoms..."
    echo "   - $native_model_pdb"
    echo "   - $native_ligand_pdb"
    echo "   - $native_ligandH_pdb"
    echo "   - $native_ligandH_mol2"
    echo ""

    for i in $(seq 1 $(grep -c ^MODEL vina_ligands.pdb)); do

        # Ligand - Convert each docked ligand pdb file: make sure it is the same as native_ligand.pdb except the coordinates
        echo "$ $python $scripts_dir/reorder_pdb_by_coord_mapping.py vina_ligand_ref.pdb vina_ligand.pdb vina_ligand_${i}.pdb > x_vina_ligand_${i}.pdb"
                $python $scripts_dir/reorder_pdb_by_coord_mapping.py vina_ligand_ref.pdb vina_ligand.pdb vina_ligand_${i}.pdb > x_vina_ligand_${i}.pdb
        echo ""
        echo "$ $python $scripts_dir/update_pdb_coord.py x_vina_ligand_${i}.pdb --ref $native_ligand_pdb -o vina_ligand_${i}_convert.pdb"
                $python $scripts_dir/update_pdb_coord.py x_vina_ligand_${i}.pdb --ref $native_ligand_pdb -o vina_ligand_${i}_convert.pdb
        echo ""
        echo "$ $python $scripts_dir/check_2pdbs.py $native_ligand_pdb vina_ligand_${i}_convert.pdb"
                $python $scripts_dir/check_2pdbs.py $native_ligand_pdb vina_ligand_${i}_convert.pdb
        echo ""
        
       ## LigandH - Convert vina_ligand_${i}_convert.pdb to vina_ligandH_${i}_convert.pdb by adding hydrogen atoms
       #echo "$ $python $scripts_dir/addh.py vina_ligand_${i}_convert.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o vina_ligandH_${i}_convert.pdb"
       #        $python $scripts_dir/addh.py vina_ligand_${i}_convert.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o vina_ligandH_${i}_convert.pdb
       #echo ""
       #echo "$ $obabel vina_ligandH_${i}_convert.pdb -O x_vina_ligandH_${i}_convert.mol2"
       #        $obabel vina_ligandH_${i}_convert.pdb -O x_vina_ligandH_${i}_convert.mol2
       #echo ""
       #echo "$ $python $scripts_dir/check_2pdbs.py $native_ligandH_pdb vina_ligandH_${i}_convert.pdb"
       #        $python $scripts_dir/check_2pdbs.py $native_ligandH_pdb vina_ligandH_${i}_convert.pdb
       #echo ""

    done

    echo "# Check their RMSD values: vina_ligand_*_convert.pdb and vina_ligandH_*_convert.pdb"
    echo "# The RMSD values: vina_ligand_*_convert.pdb"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligand_*_convert.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligand_*_convert.pdb"
    echo ""

   #echo "# The RMSD values: vina_ligandH_*_convert.pdb"
   #echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligandH_*_convert.pdb""
   #        $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligandH_*_convert.pdb"
   #echo ""

else # using the default $input_pdb

    native_model_pdb="$work_dir/$input_pdb_base"
    native_ligand_pdb="$vina_dir/vina_ligand.pdb"
    echo "# ✅ Using the input_pdb as the reference model structure. Proceeding with ligand alignment..."
    echo "   - $native_model_pdb"

    echo "# Convert each docked ligand pdb file and make sure it is the same as native_ligand.pdb except the coordinates"
    for i in $(seq 1 $(grep -c ^MODEL vina_ligands.pdb)); do

        echo "$ $python $scripts_dir/reorder_pdb_by_coord_mapping.py vina_ligand_ref.pdb x_ligand.pdb vina_ligand_${i}.pdb > x_vina_ligand_${i}.pdb"
                $python $scripts_dir/reorder_pdb_by_coord_mapping.py vina_ligand_ref.pdb x_ligand.pdb vina_ligand_${i}.pdb > x_vina_ligand_${i}.pdb
        echo ""
        echo "$ $python $scripts_dir/update_pdb_coord.py x_vina_ligand_${i}.pdb --ref $native_ligand_pdb -o vina_ligand_${i}_convert.pdb"
                $python $scripts_dir/update_pdb_coord.py x_vina_ligand_${i}.pdb --ref $native_ligand_pdb -o vina_ligand_${i}_convert.pdb
        echo ""
        echo "$ $python $scripts_dir/check_2pdbs.py $native_ligand_pdb vina_ligand_${i}_convert.pdb"
                $python $scripts_dir/check_2pdbs.py $native_ligand_pdb vina_ligand_${i}_convert.pdb
        echo ""
        
    done

    echo "# Check their RMSD values: aligned_model.pdb & aligned_ligand.pdb:"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligand_*_convert.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "vina_ligand_*_convert.pdb"
    echo ""

fi

echo "# Remove the temporary files"
echo "$ rm x_vina*.pdb"
        rm x_vina*.pdb
echo ""

echo "# Go back to the current directory"
echo "$ cd $default_work_dir"
        cd $default_work_dir
echo ""

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for vina.sh: $ELAPSED_TIME seconds"

