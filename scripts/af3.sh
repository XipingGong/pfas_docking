#!/bin/bash
#SBATCH --job-name=af3			#Name your job something original
#SBATCH --partition=gpu_p			#Use the GPU partition
#SBATCH --ntasks=1		
#SBATCH --cpus-per-task=32			#If you use the default options, AlphaFold3 will run four simutaneous Jackhmmer processes with 8 CPUs each
#SBATCH --gres=gpu:1				#If you don’t care whether your job uses an A100 node or an H100 node (and there isn’t much difference in run time)…
#SBATCH --constraint=Milan|SapphireRapids	#…this is the easiest way to specify either one without accidentally using a P100 or L4, which lack sufficient device memory
#SBATCH --mem=60gb
#SBATCH --time=120:00:00
#SBATCH --output=%x.%j.out     
#SBATCH --error=%x.%j.out  # Redirect stderr to the same file as stdout 

SECONDS=0  # Reset the counter

work_dir=`pwd` # it must have the input.pdb (required)

scripts_dir='/home/xg69107/work/pfas_pdbs/pfas_docking/scripts' # need to change
af3_param_dir='/home/xg69107/program/alphafold3' # need to change ; the directory to include the AlphaFold 3 model parameters

# unnecessary to change unless you know what you are doing
echo "$ python $scripts_dir/get_json_for_af3.py $work_dir/input.pdb"
        python $scripts_dir/get_json_for_af3.py $work_dir/input.pdb
echo ""

echo "$ cat $work_dir/input.json"
        cat $work_dir/input.json
echo ""

echo "Running af3 model >>"
singularity exec \
     --nv \
     --bind $work_dir:/root/af_input \
     --bind $work_dir:/root/af_output \
     --bind $af3_param_dir:/root/models \
     --bind /db/AlphaFold3/20241114:/root/public_databases \
     /apps/singularity-images/alphafold-3.0.0-CCDpatched.sif \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_input/input.json \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output

# Create the directory for best docking pose
echo "$ mkdir -p $work_dir/input/best_pose"
        mkdir -p $work_dir/input/best_pose
echo "$ cp $work_dir/input/input_model.cif $work_dir/input/best_pose/model.cif"
        cp $work_dir/input/input_model.cif $work_dir/input/best_pose/model.cif

# Post-processing
echo "$ python $scripts_dir/convert_cif_to_pdb.py $work_dir"
        python $scripts_dir/convert_cif_to_pdb.py $work_dir
echo ""

echo "$ python $scripts_dir/align_model_and_ligand.py --dir $work_dir $work_dir/input.pdb"
        python $scripts_dir/align_model_and_ligand.py --dir $work_dir $work_dir/input.pdb
echo ""

echo "$ python $scripts_dir/combine_aligned_ligands.py --dir $work_dir -o $work_dir/af3_ligands.pdb"
        python $scripts_dir/combine_aligned_ligands.py --dir $work_dir -o $work_dir/af3_ligands.pdb
echo ""

echo "$ python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/af3_ligands.pdb"
        python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/af3_ligands.pdb
echo ""

echo "# Elapsed time: $SECONDS seconds"

