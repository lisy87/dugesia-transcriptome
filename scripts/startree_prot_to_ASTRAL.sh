##script_name: startree_prot_to_ASTRAL.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N iqtree_ft             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate iqtree

THREADS=$1

for FILE in *.fasta
do

iqtree -s $FILE -m LG+F+G -nt $THREADS

# iqtree -s $FILE -m MFP -bb 10000 -nt $THREADS Also, you can use this command instead

echo "*****************" iqtree $FILE finished "************************************"
echo "*******************************************************************************"
done

