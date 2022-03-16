##script_name: PhyloBayes_PROT.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N phylobayes             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate cc_phylobayes

THREADS=$1
INFILE=$2 #phylip format
MODEL=$3 # -lg for prot
INTREE=$4 #starting tree
CHAINNAME=$5 # submit two jobs: chain_o and chain_1

mpirun -n $THREADS pb_mpi -d $INFILE -catfix C20 $MODEL -t $INTREE $CHAINNAME

echo PhyloBayes protein finished
