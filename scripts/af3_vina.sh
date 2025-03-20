#!/bin/bash

START_TIME=$SECONDS  # Start timer for this PDB file
scripts_dir='/home/xg69107/work/pfas_pdbs/pfas_docking/scripts' # need to change

echo "python $scripts_dir/get_inputs_for_vina.py input.pdb"
python $scripts_dir/get_inputs_for_vina.py input.pdb
echo ""
    
# Keep the parameter file, but replace the protein and ligand pdbs with native ones 
cp ../../x_receptor.pdb `pwd` # copy the native protein here, to replace the predicted one
cp ../../x_ligand.pdb `pwd` # copy the native ligand here, to replace the predicted one

echo "mk_prepare_receptor.py -i x_receptor.pdb -o x_receptor -p -a"
mk_prepare_receptor.py -i x_receptor.pdb -o x_receptor -p -a
echo ""

# Prepare Ligand PDB

obabel x_ligand.pdb -O x_ligand.sdf -p 7.4 # pH 7.4
echo ""

echo "mk_prepare_ligand.py -i x_ligand.sdf -o x_ligand.pdbqt"
mk_prepare_ligand.py -i x_ligand.sdf -o x_ligand.pdbqt
echo ""

obabel x_ligand.pdbqt -O vina_ligand_ref.pdb # as a reference
echo ""

# Run AutoDock Vina docking
vina --receptor x_receptor.pdbqt \
     --ligand x_ligand.pdbqt \
     --config x_box_params.txt \
     --out x_ligand_docked.pdbqt \
     --exhaustiveness 32 \
     --num_modes 5
echo ""

# Convert docked ligand output from PDBQT to PDB
obabel x_ligand_docked.pdbqt -O vina_ligands.pdb --separate
echo ""

python $scripts_dir/rmsd.py vina_ligand_ref.pdb vina_ligands.pdb

END_TIME=$SECONDS  # End timer for this PDB file
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "$(date '+%Y-%m-%d %H:%M:%S') Time taken for $PDB_FILE: $ELAPSED_TIME seconds"

