#!/bin/bash

SCRIPTS_DIR="/home/xg69107/work/pfas_pdbs/pfas_docking/scripts"

PDB_FILE=$1 # /xx/pdbs/6WJ5_LXY.pdb
PDB_DIR=$2  # /xx/dock_dir/6WJ5_LXY

PDB_FILENAME=$(basename "$PDB_FILE")
[[ "$PDB_FILENAME" != *.pdb ]] && { echo "Error: Input file must be a .pdb file."; exit 1; }

IFS='_' read -r PDBID LIGAND_ID <<< "${PDB_FILENAME%.pdb}"
[ -z "$PDBID" ] || [ -z "$LIGAND_ID" ] && { echo "Error: Filename format should be <PDBID>_<LIGAND>.pdb"; exit 1; }

echo "#---------------------"
echo "$ mkdir -p "$PDB_DIR""
        mkdir -p "$PDB_DIR"

echo "#---------------------"
echo "$ cp $PDB_FILE $PDB_DIR"
        cp $PDB_FILE $PDB_DIR

echo "#---------------------"
echo "$ cd $PDB_DIR"
        cd $PDB_DIR

PDB_NAME=${PDBID}_${LIGAND_ID}
echo "#---------------------"
echo "$ python "$SCRIPTS_DIR/check_pdb.py" "${PDB_NAME}.pdb" --ligand_id "$LIGAND_ID""
        python "$SCRIPTS_DIR/check_pdb.py" "${PDB_NAME}.pdb" --ligand_id "$LIGAND_ID"

echo "#---------------------"
echo "$ python "$SCRIPTS_DIR/clean_pdb.py" "${PDB_NAME}.pdb" -o "${PDB_NAME}_cleaned.pdb""
        python "$SCRIPTS_DIR/clean_pdb.py" "${PDB_NAME}.pdb" -o "${PDB_NAME}_cleaned.pdb"

echo "#---------------------"
echo "$ python "$SCRIPTS_DIR/extract_pdb.py" "${PDB_NAME}_cleaned.pdb" "$LIGAND_ID""
        python "$SCRIPTS_DIR/extract_pdb.py" "${PDB_NAME}_cleaned.pdb" "$LIGAND_ID"

echo "#---------------------"
echo "$ mv *Protein*.pdb input.pdb"
        mv *Protein*.pdb input.pdb

echo "#---------------------"
echo "$ python "$SCRIPTS_DIR/check_pdb.py" "input.pdb" --ligand_id "$LIGAND_ID""
        python "$SCRIPTS_DIR/check_pdb.py" "input.pdb" --ligand_id "$LIGAND_ID"

echo "#---------------------"
echo "$ ls -lrt `pwd`"
        ls -lrt `pwd`

echo "#---------------------"
echo "Processing completed for $PDB_NAME"



