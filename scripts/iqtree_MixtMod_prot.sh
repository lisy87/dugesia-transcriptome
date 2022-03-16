##script_name: iqtree_MixtMod_prot.sh

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
MODEL=$2 # LG+C20+F+G 
         # You can use categories C20-C60 or use the PMSF aproximation
TREE=$3 # starting tree
REP=$4 #-bb Ultrafast Bootstrap I recomend a minimun of 1000000 replicates
THREADS=$5

iqtree -s $FILE -m $MODEL -ft $TREE -bb $REP -nt $THREADS

echo iqtree finished

#NOTE: with LG+C20+F+G, -bb 1000000, 162232 positions, and 83 samples:
#CPU time used for tree search: 772525.416 sec (214h:35m:25s)
#Wall-clock time used for tree search: 51815.150 sec (14h:23m:35s)
#Total CPU time used: 790585.844 sec (219h:36m:25s)
#Total wall-clock time used: 62278.360 sec (17h:17m:58s)

