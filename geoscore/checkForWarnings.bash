#!/bin/bash

# Check warnings in code (and seperate ours from those in external libs)

echo 'Cleaning'
##############

rm -f make_output_*.txt warnings_all_*.txt warnings_geos_*.txt

make clean

echo 'Making (GCC)'
##################

./make_gcc.bash &> make_output_gcc.txt

echo 'Getting warnings'
#######################

grep warning make_output_gcc.txt > warnings_all_gcc.txt

grep -v 'external_libs' warnings_all_gcc.txt  > warnings_geos_gcc.txt


echo 'Cleaning'
##############

make clean

echo 'Making (ICC)'
##################

./make_icc.bash &> make_output_icc.txt

echo 'Getting warnings'
#######################

grep warning make_output_icc.txt > warnings_all_icc.txt

grep -v 'external_libs' warnings_all_icc.txt  > warnings_geos_icc.txt

echo 'Done'
###########
