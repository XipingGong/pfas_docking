#!/bin/bash
#SBATCH --job-name=af3 #Name your job something original
#SBATCH --partition=gpu_p #Use the GPU partition
##SBATCH --partition=gpu_30d_p #Use the 30-day GPU partition
#SBATCH --ntasks=1	
#SBATCH --cpus-per-task=32 #If you use the default options, AlphaFold3 will run four simutaneous Jackhmmer processes with 8 CPUs each
#SBATCH --gres=gpu:1 #If you don’t care whether your job uses an A100 node or an H100 node (and there isn’t much difference in run time)…
#SBATCH --constraint=Milan|SapphireRapids #…this is the easiest way to specify either one without accidentally using a P100 or L4, which lack sufficient device memory
#SBATCH --mem=60gb
#SBATCH --output=%x.%j.out     
#SBATCH --time=0-01:00:00 # request 1-hour to run the AF3 job
#SBATCH --error=%x.%j.out  # Redirect stderr to the same file as stdout 

START_TIME=$SECONDS

# This bash script "af3.sh" is used to run AlphaFold 3 for protein-ligand docking.
# This script has been successfully tested in the computing center @UGA.
# To run it, you may need to change the INPUT section in this script.

# INPUT
# -----

# This is a "scripts" folder that can be downloaded from https://github.com/XipingGong/pfas_docking.git
scripts_dir='/home/xg69107/program/pfas_docking/scripts' # need to change 

# This is a folder that has the AF3 parameters file requested from DeepMind
# The parameters file needs to be requested, and please see the link: https://github.com/google-deepmind/alphafold3
af3_param_dir='/home/xg69107/program/alphafold3' # need to change

# These should be installed
obabel="/home/xg69107/program/anaconda/anaconda3/envs/gmxMMPBSA/bin/obabel"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"

# INFO
# ----
print_help() {
    echo ""
    echo "Usage: sbatch af3.sh --input_json FILE [--work_dir DIR] [--native_dir DIR]"
    echo ""
    echo "Required:"
    echo "  --input_json     JSON file containing protein (chain A) and ligand (chain B)"
    echo ""
    echo "Optional:"
    echo "  --work_dir       Working directory (default: current directory)"
    echo "  --native_dir     Folder that contains the reference structures, such as native_model.pdb, native_ligand.pdb, native_ligandH.pdb, and native_ligandH.mol2,"
    echo "                   (default: if unprovided, then the predicted best pose will be a reference for the RMSD calculations)"
    echo "  --run_af3        Whether to run the AF3 docking command (default: true). Set to false to skip the docking step."
    echo "  -h, --help       Show this help message and exit"
    echo ""
    exit 0
}

# Default Values 
default_work_dir=$(pwd)
work_dir=$default_work_dir
native_dir="$work_dir/af3/best_pose"
run_af3=true  # default is true

# This JOSN file has the info of one protein and one ligand, and their residue name must be in the RCSB PDB.
# This file must be in the working directory.
# It must have two chains, like chain A for protein and chain B for ligand.
input_json=""

# Parse Arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_json) input_json="$2"; shift ;;
        --work_dir) work_dir="$2"; shift ;;
        --native_dir) native_dir="$2"; shift ;;
        --run_af3) run_af3="$2"; shift ;;
        -h|--help) print_help ;;
        *) echo "❌ Unknown parameter: $1"; print_help ;;
    esac
    shift
done
if [ -z "$input_json" ]; then
    echo "❌ Error: --input_json is required."
    print_help
fi

# Convert input_json to full path
[[ "$input_json" != /* ]] && input_json="$default_work_dir/$input_json"

# Validate $input_json
if [ ! -f "$input_json" ]; then
    echo "❌ Error: $input_json must be a valid file."
    print_help
fi
input_json_base="${input_json##*/}" # e.g., af3.json
input_json_basename="${input_json_base%.*}" # e.g., af3

# Convert work_dir to full path
[[ "$work_dir" != /* ]] && work_dir="$default_work_dir/$work_dir"

# Validate $work_dir
if [ ! -d "$work_dir" ]; then
    echo "❌ Error: work_dir does not exist: $work_dir"
    exit 1
fi

# Convert native_dir to full path
[[ "$native_dir" != /* ]] && native_dir="$default_work_dir/$native_dir"

# OPTIONAL (you could need to modify AF3 section, but do not modify others, except you know what you are doing)
# --------

# Starting
echo "# Running an AF3 job >>"
echo ""

# Go to the $work_dir directory
echo "$ cd $work_dir"
        cd $work_dir
echo ""

# Create the af3_dir for saving the data
af3_dir="$work_dir/$input_json_basename"

# Validate $af3_dir
if [ -d "$af3_dir" ] && [ "$(ls -A "$af3_dir")" ]; then
    echo "⚠️  Warning: AF3 folder already exists and is not empty:"
    echo "   → $af3_dir"
    echo "   Existing data may be overwritten."
else
    echo "✅ All data will be saved into this AF3 folder: $af3_dir"
    echo "mkdir -p $af3_dir"
          mkdir -p $af3_dir
    echo ""
fi

# Copy $input_json in the $work_dir, and show the $input_json file for AF3
if [ -f "$input_json_base" ]; then
    if ! cmp -s "$input_json" "$input_json_base"; then
        echo "❌ A different $input_json_base already exists in the current directory!"
        echo "   Existing: $(readlink -f "$input_json_base")"
        echo "   Provided: $(readlink -f "$input_json")"
        exit 1
    else
        echo "ℹ️  $input_json_base already exists in the $work_dir and is identical. Skipping copy."
        echo ""
    fi
else
    echo "$ cp $input_json $input_json_base"
            cp $input_json $input_json_base
    echo ""
fi

echo "# AF3 input: $input_json_base"
echo "$ cat $input_json_base"
        cat $input_json_base
echo ""

# Run the AF3 docking (change it if needed)
if [[ "$run_af3" == "true" ]]; then
    echo "# Run the AF3 docking"
    echo "$ singularity exec --nv --bind $work_dir:/root/af_input --bind $work_dir:/root/af_output --bind $af3_param_dir:/root/models --bind /db/AlphaFold3/20241114:/root/public_databases /apps/singularity-images/alphafold-3.0.0-CCDpatched.sif python /app/alphafold/run_alphafold.py --json_path=/root/af_input/$input_json_base --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output"
    singularity exec \
        --nv \
        --bind $work_dir:/root/af_input \
        --bind $work_dir:/root/af_output \
        --bind $af3_param_dir:/root/models \
        --bind /db/AlphaFold3/20241114:/root/public_databases \
        /apps/singularity-images/alphafold-3.0.0-CCDpatched.sif \
        python /app/alphafold/run_alphafold.py \
        --json_path=/root/af_input/$input_json_base \
        --model_dir=/root/models \
        --db_dir=/root/public_databases \
        --output_dir=/root/af_output
else
    echo "# ⏩ Skipping AF3 docking step (singularity exec) as --run_af3 is set to false."
fi
echo ""

# Analysis
# --------

echo "# Create a "best_pose" folder & copy the af3_model.cif file into model.cif in this folder"
echo "$ mkdir -p $af3_dir/best_pose"
        mkdir -p $af3_dir/best_pose
echo "$ cp $af3_dir/${input_json_basename}_model.cif $af3_dir/best_pose/model.cif"
        cp $af3_dir/${input_json_basename}_model.cif $af3_dir/best_pose/model.cif
echo ""

# Convert cif to pdb
for folder in "$af3_dir"/*; do
    [ -f "$folder/model.cif" ] || continue
    echo "$ $obabel $folder/model.cif -O $folder/model.pdb"
            $obabel $folder/model.cif -O $folder/model.pdb
done
echo ""

echo "# Align the structures with the native reference structures"
echo "# Why this step: the predicted structure files could miss the atom info included in the reference files."
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

    echo "# Convert models - Convert model.pdb to model_convert.pdb and make sure it is the same format as the reference except the coordinates"
    for folder in "$af3_dir"/[bs]*; do
        echo "$ $python $scripts_dir/reorder_pdb_by_atom_info.py --ref $native_model_pdb $folder/model.pdb -o $folder/model_convert.pdb"
                $python $scripts_dir/reorder_pdb_by_atom_info.py --ref $native_model_pdb $folder/model.pdb -o $folder/model_convert.pdb
        echo ""
        echo "$ $python $scripts_dir/check_2pdbs.py $native_model_pdb $folder/model_convert.pdb"
                $python $scripts_dir/check_2pdbs.py $native_model_pdb $folder/model_convert.pdb
    done
    echo ""

    echo "# Alignment - Create the aligned the predicted model & ligand PDB structures"
    echo "$ $python $scripts_dir/align_pdb.py "$af3_dir/[bs]*/model_convert.pdb" --ref "$native_model_pdb""
            $python $scripts_dir/align_pdb.py "$af3_dir/[bs]*/model_convert.pdb" --ref "$native_model_pdb"
    echo ""

    echo "# Rename algned_model.pdb & aligned_ligand.pdb"
    for folder in "$af3_dir"/[bs]*; do
        echo "mv $folder/aligned_model.pdb $folder/aligned_model_convert.pdb"
              mv $folder/aligned_model.pdb $folder/aligned_model_convert.pdb
        echo "mv $folder/aligned_ligand.pdb $folder/aligned_ligand_convert.pdb"
              mv $folder/aligned_ligand.pdb $folder/aligned_ligand_convert.pdb
    done
    echo ""

   #echo "# LigandH - Convert aligned_ligand_convert.pdb to aligned_ligandH_convert.pdb by adding hydrogen atoms"
   #for folder in "$af3_dir"/[bs]*; do
   #    echo "$ $python $scripts_dir/addh.py $folder/aligned_ligand_convert.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o $folder/aligned_ligandH_convert.pdb"
   #            $python $scripts_dir/addh.py $folder/aligned_ligand_convert.pdb --ref_pdb $native_ligandH_pdb --ref_mol2 $native_ligandH_mol2 -o $folder/aligned_ligandH_convert.pdb
   #    echo ""
   #    echo "$ $obabel $folder/aligned_ligandH_convert.pdb -O $folder/x_aligned_ligandH_convert.mol2"
   #            $obabel $folder/aligned_ligandH_convert.pdb -O $folder/x_aligned_ligandH_convert.mol2
   #    echo ""
   #    echo "$ $python $scripts_dir/check_2pdbs.py $native_ligandH_pdb $folder/aligned_ligandH_convert.pdb"
   #            $python $scripts_dir/check_2pdbs.py $native_ligandH_pdb $folder/aligned_ligandH_convert.pdb
   #    echo ""
   #done

    echo "# Check the RMSD values"
    echo "# The RMSD values: aligned_model_convert.pdb"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_model_convert.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_model_convert.pdb"
    echo ""

    echo "# The RMSD values: aligned_ligand_convert.pdb"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligand_convert.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligand_convert.pdb"
    echo ""

   #echo "# The RMSD values: aligned_ligandH_convert.pdb"
   #echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligandH_convert.pdb""
   #        $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligandH_convert.pdb"
   #echo ""

else # using the default AF3-predicted best pose

    echo "# Alignment - Create the aligned the predicted model & ligand PDB structures"
    native_model_pdb="$af3_dir/best_pose/model.pdb"
    echo "# ✅ Using the AF3-predicted best pose as the reference model structure. Proceeding with ligand alignment..."
    echo "   - $native_model_pdb"
    echo "$ $python $scripts_dir/align_pdb.py "$af3_dir/[bs]*/model.pdb" --ref "$native_model_pdb""
            $python $scripts_dir/align_pdb.py "$af3_dir/[bs]*/model.pdb" --ref "$native_model_pdb"
    echo ""

    echo "# Check their RMSD values: aligned_model.pdb"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_model.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_model.pdb"
    echo ""

    echo "# Check their RMSD values: aligned_ligand.pdb"
    echo "$ $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligand.pdb""
            $python $scripts_dir/check_rmsd.py --ref $native_model_pdb "$af3_dir/[bs]*/aligned_ligand.pdb"
    echo ""

    # add hydrogens

fi

echo "# Remove the temporary files"
echo "$ rm -f $af3_dir/[bs]*/x_*.pdb"
        rm -f $af3_dir/[bs]*/x_*.pdb
echo ""

echo "# Go back to the default work_dir "
echo "$ cd $default_work_dir"
        cd $default_work_dir
echo ""

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for af3.sh: $ELAPSED_TIME seconds"


