#!/bin/bash

# Figure 1
echo ">> Figure 1"
echo "-------------"
cp before_set.dat x.log
python ../scripts/rmsd_stat.py 'x.log' --columns '8'  --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '12' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '14' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '16' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '18' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '20' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '22' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '24' --range '[0,0.2)'

cp after_set.dat x.log
python ../scripts/rmsd_stat.py 'x.log' --columns '8'  --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '12' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '14' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '16' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '18' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '20' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '22' --range '[0,0.2)'
python ../scripts/rmsd_stat.py 'x.log' --columns '24' --range '[0,0.2)'

echo ""
# exit

# Figure 2

echo ">> Figure 2"
echo "-------------"
cp before_set.dat x.log
awk '{if ($8 >= 0.2) print}' x.log > x1.log
python ../scripts/rmsd_stat.py x1.log --columns '9 and 11' --range '[0,0.2)'
python ../scripts/rmsd_stat.py x1.log --columns '9'       --range '[0.2,1000)'
python ../scripts/rmsd_stat.py x1.log --columns '11'       --range '[0.2,1000)'
python ../scripts/rmsd_stat.py x1.log --columns '9 and 11' --range '[0.2,1000)'

cp after_set.dat x.log
awk '{if ($8 >= 0.2) print}' x.log > x1.log
python ../scripts/rmsd_stat.py x1.log --columns '9 and 11' --range '[0,0.2)'
python ../scripts/rmsd_stat.py x1.log --columns '9'       --range '[0.2,1000)'
python ../scripts/rmsd_stat.py x1.log --columns '11'       --range '[0.2,1000)'
python ../scripts/rmsd_stat.py x1.log --columns '9 and 11' --range '[0.2,1000)'

echo ""
#exit

# Figure 3
echo ">> Figure 3"
echo "-------------"
awk '{print $1, $3, $4, $5}' before_set.dat > time_before.dat
awk '{print $1, $3, $4, $5}' after_set.dat > time_after.dat
echo "Data are saved into time_before.dat and time_after.dat"

echo ""
#exit


# Figure 4
echo ">> Figure 4"
echo "-------------"
python ../scripts/rmsd_hist.py --columns '8 14' before_set.dat > rmsd_hist_before.dat
python ../scripts/rmsd_hist.py --columns '8 14' after_set.dat > rmsd_hist_after.dat
echo "Data are saved into the rmsd_hist_before.dat and rmsd_hist_after.dat"

echo ""
#exit

# Visualization
echo ">> Generate data for visualization, to identify the PDBID list for Table 1 and 2"
echo "-------------"
awk '{split($1, a, "_"); print a[2]}' before_set.dat > x.log
python ../scripts/extract_ligand_info.py ccdcif_info.log x.log > visual_before.dat

awk '{split($1, a, "_"); print a[2]}' after_set.dat > x.log
python ../scripts/extract_ligand_info.py ccdcif_info.log x.log > visual_after.dat

echo "Data are saved into the visual_before.dat and visual_after.dat"
echo ""

#exit

# Table 1
echo ">> Table 1"
echo "-------------"

# R-CF2-CF(R')(R") where R, R', R" ≠ H; US EPA OPPT definition
smile='[CX4]([!#1])(F)(F)-[CX4](F)([!#1])([!#1])'
python ../scripts/filter_ligands_by_smile.py visual_after.dat $smile 
python ../scripts/filter_ligands_by_smile.py visual_before.dat $smile

echo ""
# - Using visual_before.dat and visual_after.dat to visualize the PFAS
# - Before Set
grep 3RZ7_RZ7 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 6RZX_KPQ before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 3RZ1_RZ1 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 7JYM_Z8I before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 6VQF_R7V before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 7JTM_VK7 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 7LUK_YDY before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 3RYZ_RYZ before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 5DDF_5A1 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 4J03_FVS before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 3RYX_RYX before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
# - After Set
grep 7FEU_4I6 after_set.dat  | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'

echo ""
#exit


# Table 2
echo ">> Table 2"
echo "-------------"
# phenolic fluorine (-Ph-F) ligands
smile='Oc1cccc(F)c1'
python ../scripts/filter_ligands_by_smile.py visual_before.dat $smile
python ../scripts/filter_ligands_by_smile.py visual_after.dat $smile 
echo ""

# - Before Set
grep B67 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 8RH before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep 42Q before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep B5R before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep BXD before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep UKB before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
grep ZI5 before_set.dat | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'
# - After Set
grep ZIS after_set.dat  | awk '{print $1, $2, $3, $8, $12, $14, $16, $18, $20, $22, $24}'

echo ""
exit

