#!/bin/bash

SECONDS=0  # Reset the counter

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

    # PDBID and Ligand ID
    PDBID=$(echo "$sub_dir" | awk -F'_' '{print $1}')
    Ligand_ID=$(echo "$sub_dir" | awk -F'_' '{print $2}')
    echo "PDBID Name: $PDBID"
    echo "Ligand_ID Name: $Ligand_ID"

    work_dir=$job_dir

    ## Run the script and save output
    echo "$ cd $job_dir"
            cd $job_dir

    # Released date
    echo "<< Extract the released date of $PDBID >>"
    echo "python $scripts_dir/get_date_for_pdbid.py $PDBID"
          python $scripts_dir/get_date_for_pdbid.py $PDBID

    # Info from af3.out
    echo "<< Extract INFO from af3.out >>"
    grep "\"sequence\":" $work_dir/af3.out
    grep "# Elapsed time" $work_dir/af3.out
    grep "Processing chain A took" $work_dir/af3.out
    grep 'mdtraj.Trajectory' $work_dir/af3.out | head -1

    # Calculate the RMSD values
    echo "<< AF3 >>"
    echo "$ python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/input/best_pose/aligned_model.pdb"
            python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/input/best_pose/aligned_model.pdb
    echo "$ python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/af3_ligands.pdb"
            python $scripts_dir/rmsd.py $work_dir/input.pdb $work_dir/af3_ligands.pdb

    echo "<< Vina >>"
    echo "$ python $scripts_dir/rmsd.py $work_dir/vina_ligand_ref.pdb $work_dir/vina_ligands.pdb"
            python $scripts_dir/rmsd.py $work_dir/vina_ligand_ref.pdb $work_dir/vina_ligands.pdb

    echo "<< AF3Pocket-Vina >>"
    echo "$ python $scripts_dir/rmsd.py $work_dir/vina_ligand_ref.pdb $work_dir/input/best_pose/vina_ligands.pdb"
            python $scripts_dir/rmsd.py $work_dir/vina_ligand_ref.pdb $work_dir/input/best_pose/vina_ligands.pdb

    echo ""

done < "$job_list"

echo "# Elapsed time: $SECONDS seconds"

