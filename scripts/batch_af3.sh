#!/bin/bash
#SBATCH --job-name=af3	#Name your job something original
#SBATCH --partition=gpu_p #Use the GPU partition
#SBATCH --ntasks=1		
#SBATCH --cpus-per-task=32 #If you use the default options, AlphaFold3 will run four simutaneous Jackhmmer processes with 8 CPUs each
#SBATCH --gres=gpu:1	#If you don’t care whether your job uses an A100 node or an H100 node (and there isn’t much difference in run time)…
#SBATCH --constraint=Milan|SapphireRapids #…this is the easiest way to specify either one without accidentally using a P100 or L4, which lack sufficient device memory
#SBATCH --mem=60gb
#SBATCH --time=120:00:00
#SBATCH --output=%x.%j.out     
#SBATCH --error=%x.%j.out  # Redirect stderr to the same file as stdout 

scripts_dir="/home/xg69107/work/pfas_pdbs/pfas_docking/scripts" # need to change
wdir='/home/xg69107/work/pfas_pdbs/pfas_docking/test' # need to change; a working directory
cd $wdir

# Read job directories from file
job_list="jobs.log"  # The file containing job directories

# Check if file exists
if [[ ! -f "$job_list" ]]; then
    echo "Error: Job directory list file '$job_list' not found!"
    exit 1
fi

# Read the job directories from the file line by line
while IFS= read -r pdb_file; do

    # Ensure the line is not empty and is a valid .pdb file
    if [[ ! "$pdb_file" =~ \.pdb$ ]]; then
        echo "Skipping invalid entry: '$pdb_file' (Not a .pdb file)"
        continue
    fi

    # Extract the base filename (e.g., "6WJ5_LXY" from "pdbs/6WJ5_LXY.pdb")
    sub_dir=$(basename "$pdb_file" .pdb)
    job_dir="$wdir/dock_dir/$sub_dir"

    echo "Processing: $sub_dir (Job Directory: $job_dir)"
    echo "-----------------------------------------"

   ## Remove all except input pdb files
   #echo "$ find dock_dir/$sub_dir ! -name "input.pdb" ! -name "$sub_dir.pdb" ! -name "${sub_dir}_cleaned.pdb" ! -name "$sub_dir" -exec rm -rf {} +"
   #        find dock_dir/$sub_dir ! -name "input.pdb" ! -name "$sub_dir.pdb" ! -name "${sub_dir}_cleaned.pdb" ! -name "$sub_dir" -exec rm -rf {} +
   #        continue

    # Create the directory and input files
    mkdir -p "$job_dir"
    echo "$ bash $scripts_dir/mk_dock_dir.sh $wdir/pdbs/${sub_dir}.pdb $job_dir > ${job_dir}/mk_dock_dir.out 2>&1"
            bash $scripts_dir/mk_dock_dir.sh $wdir/pdbs/${sub_dir}.pdb $job_dir > ${job_dir}/mk_dock_dir.out 2>&1

    ## Go to the created directory
    echo "$ cd $job_dir"
            cd $job_dir

    # Copy the af3.sh script
    echo "$ cp $scripts_dir/af3.sh $job_dir"
            cp "$scripts_dir/af3.sh" "$job_dir"

    echo "$ bash af3.sh > af3.out 2>&1"
            bash af3.sh > af3.out 2>&1

    ## Display the results, e.g., the RMSD values
    echo "$ grep 'RMSD' af3.out"
            grep 'RMSD' af3.out
    
    grep 'Elapsed time' af3.out

    echo ""

done < "$job_list"

