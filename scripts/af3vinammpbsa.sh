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

# 
# INFO
# ----
# A general procedure is as follows,
# 
# **Step 1**: Run the AF3 from the sequence information, and collect the top-5 poses
# 
# **Step 2**: Take the best pose predicted by AF3 to run the Vina docking, and then collect the top-5 poses
# 
# **Step 3**: Use the "mmpbsa.sh" script to obtain their binding free energies for the poses predicted

#
# INPUT
# -----
scripts_dir="/home/xg69107/program/pfas_docking/scripts"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"
work_dir="/home/xg69107/program/pfas_docking/test/dock_dir/7AAI_8PF" # the current working directory

af3_json="$work_dir/af3.json" # we need an "af3.json" file in the working directory
if [[ ! -f "$af3_json" ]]; then
    echo "❌ Error: Required file '$af3_json' not found. Exiting."
    exit 1
fi

#
# You do not need to modify below, except you know what you are doing
cd $work_dir
#
# AF3
# ----
mkdir -p af3 # create an "af3" folder
bash $scripts_dir/af3.sh --input_json af3.json # it will run an AF3 job
ls -lrt af3/*/aligned_model_convert.pdb # Check out the predicted structures
# >> af3/best_pose/aligned_model_convert.pdb
# >> af3/seed-1_sample-0/aligned_model_convert.pdb
# >> af3/seed-1_sample-1/aligned_model_convert.pdb
# >> af3/seed-1_sample-2/aligned_model_convert.pdb
# >> af3/seed-1_sample-3/aligned_model_convert.pdb
# >> af3/seed-1_sample-4/aligned_model_convert.pdb
echo ""
#
# Vina
# ----
mkdir -p vina
bash $scripts_dir/vina.sh --input_pdb af3/best_pose/aligned_model_convert.pdb --work_dir vina # run a Vina job
for ligand_pdb in $(ls vina/vina_ligand_*_convert.pdb | sort); do
    num=$(basename "$ligand_pdb" | sed -E 's/.*_ligand_([0-9]+)_convert\.pdb/\1/')
    cat vina/vina_receptor.pdb "$ligand_pdb" | grep 'ATOM' > vina/vina_model_${num}_convert.pdb
done
ls -lrt vina/vina_model_*_convert.pdb
# >> vina/vina_model_1_convert.pdb
# >> vina/vina_model_2_convert.pdb
# >> vina/vina_model_3_convert.pdb
# >> vina/vina_model_4_convert.pdb
# >> vina/vina_model_5_convert.pdb
echo ""
# 
# MMPBSA
# ------
# We can see that both AF3 and Vina predict 5 condidate poses, but these predicted structures do not have the explicit hydrogen atoms. We should add the hydrogen atoms before running the MMPBSA.
# To add the hydrogen atoms, we can first pick up one pose as the native structure, for example, "af3/best_pose/aligned_model_convert.pdb", to create the native structural information. 
# However, if we can obtain the native structure from the RCSB PDB, then we should use this structure as the native structure, and we can create the native structural information before running AF3 and Vina, so that we can use this information for both AF3 and Vina by using "--native_dir [DIR]"
#
# Creat the native folder >>
mkdir -p native
cp af3/best_pose/aligned_model_convert.pdb native/native_model.pdb # pick up the best pose
bash $scripts_dir/get_ref_for_af3vinammpbsa.sh --input_pdb native/native_model.pdb --work_dir native # create the native structures and topology information
# 
# Loop all predicted structures by both AF3 and Vina
mkdir -p mmpbsa
echo "🚀 Starting MMPBSA calculations..."
for pdb_file in {native/native_model.pdb,af3/*/aligned_model_convert.pdb,vina/vina_model_*_convert.pdb}; do
    if [[ -f "$pdb_file" ]]; then

        echo ">> Running: $pdb_file"

        pdb_base_file=$(echo "$pdb_file" | sed 's|/|-|g')
        pdb_base_name="${pdb_base_file%.pdb}"

        cp $pdb_file "mmpbsa/$pdb_base_file"

        # Step 1: Generate H-added input pdb file --> 
        bash $scripts_dir/get_ref_for_af3vinammpbsa.sh \
            --input_pdb "mmpbsa/${pdb_base_name}.pdb" \
            --work_dir mmpbsa \
            --native_dir native # --> generate "mmpbsa/${pdb_base_name}H.pdb"
        
        # Step 2: Run MMPBSA
        bash $scripts_dir/mmpbsa.sh \
            --input_pdb "mmpbsa/${pdb_base_name}H.pdb" \
            --work_dir mmpbsa \
            --native_dir native

        # Step 3: Check RMSD values
        if [[ $pdb_base_name == af3* ]]; then
            ref_file="mmpbsa/af3-best_pose-aligned_model_convertH_emin.pdb"
        elif [[ $pdb_base_name == vina* ]]; then
            ref_file="mmpbsa/vina-vina_model_1_convertH_emin.pdb"
        elif [[ $base_name == nat* ]]; then
            ref_file="native/native_modelH.pdb"
        else
            echo "⚠️ Warning: Unknown prefix for $pdb_base_name — skipping RMSD"
            return
        fi
        $python $scripts_dir/check_rmsd.py \
            --ref "$ref_file" \
            "mmpbsa/${pdb_base_name}H_emin.pdb" \
            > "mmpbsa/${pdb_base_name}H_emin_Top1RMSD.dat"
        echo ""
        
    fi
        
done
#
# A quick analysis
# ----------------
# Output MMPBSA & RMSD results
echo "grep 'ΔTOTAL' mmpbsa/[nav][afi]*_MMPBSA.dat"
      grep 'ΔTOTAL' mmpbsa/[nav][afi]*_MMPBSA.dat
echo ""

echo "grep 'Ligand RMSD (Direct)' mmpbsa/[nav][afi]*_Top1RMSD.dat"
      grep 'Ligand RMSD (Direct)' mmpbsa/[nav][afi]*_Top1RMSD.dat
echo ""
#
# Identify the best pose 
grep "ΔTOTAL\|Ligand RMSD (Direct)" mmpbsa/[av][fi]*.dat > mmpbsa/run.dat
$python $scripts_dir/af3vinammpbsa_ana.py mmpbsa/run.dat

END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for af3vinammpbsa.sh script: $ELAPSED_TIME seconds"

