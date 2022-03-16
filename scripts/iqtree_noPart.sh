##script_name: iqtree_noPart.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N iqtree_fast             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate iqtree

FILE=$1
BOOTS=$2 # you can use 10000 for a first aproximation
THREADS=$3

iqtree -s $FILE -m MFP -bb $BOOTS -nt $THREADS

echo iqtree without partitions finished

