#!/bin/bash
#SBATCH --job-name=af3_vina                   # Job name
#SBATCH --partition=batch               # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=4               # Run on a single CPU
#SBATCH --mem=4GB                       # Job memory request
#SBATCH --time=3-00:00:00              # day-hour:min:seconds
#SBATCH --output=%x.%j.out         # Standard output log
#SBATCH --error=%x.%j.out          # Standard error log

scripts_dir="/home/xg69107/work/pfas_pdbs/pfas_docking/scripts" # need to change
wdir='/home/xg69107/work/pfas_pdbs/pfas_docking/test' # need to change; a working directory
cd $wdir

# Read job directories from file
job_list="jobs.log"  # The file containing job directories

# Unnecessary to change the following except you know what you are doing
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
    job_dir="$wdir/dock_dir/$sub_dir/input/best_pose"

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

    # Create the input.pdb
    if [ -f aligned_model.pdb ]; then
        echo "$ cp aligned_model.pdb input.pdb"
                cp aligned_model.pdb input.pdb
    else
        echo "${sub_dir}: File aligned_model.pdb not found. Skipping."
        echo ""
        continue
    fi

    # Copy the new script
    echo "$ cp $scripts_dir/af3_vina.sh $job_dir"
            cp $scripts_dir/af3_vina.sh $job_dir

    echo "$ bash af3_vina.sh > af3_vina.out 2>&1"
            bash af3_vina.sh > af3_vina.out 2>&1

    ## Display the results, e.g., RMSD values
    echo "$ grep 'RMSD' af3_vina.out"
            grep 'RMSD' af3_vina.out
    
    grep 'Elapsed time' af3_vina.out

    echo ""

done < "$job_list"

