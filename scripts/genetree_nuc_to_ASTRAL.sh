##script_name: genetree_nuc_to_ASTRAL.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N iqtree_ft             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate iq3

THREADS=$1

for FILE in *.fasta
do

iqtree -s $FILE -m "MIX{JC,HKY,GTR}+G4" -bb 10000 -nt $THREADS

echo "*****************" iqtree $FILE finished "************************************"
echo "*******************************************************************************"
done
