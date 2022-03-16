##script_name: PhyloBayes_NUC.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N phylobayes_nuc             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate cc_phylobayes

THREADS=$1
INFILE=$2 #phylip format
INTREE=$3 #starting tree
CHAINNAME=$4 # submit two jobs: chain_o and chain_1

  #MODEL = -cat -gtr for DNA

mpirun -n $THREADS pb_mpi -d $INFILE -cat -gtr -t $INTREE $CHAINNAME

echo PhyloBayes nucleotide finished
