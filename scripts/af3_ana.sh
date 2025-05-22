
scripts_dir="/home/xg69107/program/pfas_docking/scripts"
python="/home/xg69107/program/anaconda/anaconda3/bin/python"

NATIVE_MODEL="native_model.pdb"
NATIVE_LIGAND="native_ligand.pdb"
awk 'substr($0,22,1)=="B"' $NATIVE_MODEL > $NATIVE_LIGAND

for bp in af3/*; do
  [ -d "$bp" ] || continue

  # split out chain A (protein) and chain B (ligand)
  awk 'substr($0,22,1)=="A"' "$bp/model.pdb" > "$bp/x_prot.pdb"
  awk 'substr($0,22,1)=="B"' "$bp/model.pdb" > "$bp/x_lig.pdb"

  echo "python $scripts_dir/rename_atom_info_by_topology.py --ref $NATIVE_LIGAND $bp/x_lig.pdb > $bp/x_lig_convert.pdb"
        python $scripts_dir/rename_atom_info_by_topology.py --ref $NATIVE_LIGAND $bp/x_lig.pdb > $bp/x_lig_convert.pdb
  echo ""

  # stitch protein + converted ligand into one PDB
  cat "$bp/x_prot.pdb" "$bp/x_lig_convert.pdb" > "$bp/model_convert.pdb"

  echo "$bp: Converted $bp/model_convert.pdb"
done
echo ""

# Alignment: default 'aligned_model.pdb' will be created for each AF3-predicted structure 
$python $scripts_dir/align_pdb.py  --ref $NATIVE_MODEL "af3/*/model_convert.pdb" 

# Alignment: Convert all AF3-predicted structures into a pdb file 'af3_model.pdb' 
$python $scripts_dir/align_pdb.py  --ref $NATIVE_MODEL "af3/*/aligned_model.pdb" -o af3_model.pdb

# Calculate RMSD in terms of native structure
$python $scripts_dir/check_rmsd.py --ref $NATIVE_MODEL "af3/*/aligned_model.pdb"

