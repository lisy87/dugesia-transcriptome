##script_name: iqtree_MixtMod_nuc.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N iqtree_mix             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate iqtree

FILE=$1 # alignment in phylip format
MODEL=$2 # "MIX{JC,HKY,GTR}+G4" for nucleotide data
REP=$3 #-bb Ultrafast Bootstrap I recomend a minimun of 1000000 replicates
THREADS=$4

iqtree -s $FILE -m $MODEL -bb $REP -nt $THREADS
