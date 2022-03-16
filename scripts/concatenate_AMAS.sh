#script_name: concatenate_AMAS.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N concat                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate

OUTconc=$1
OUTpart=$2
TYPE=$3 #(dna,aa)
CORES=$4

AMAS.py concat -i *.fasta -f fasta -d $TYPE -p $OUTpart -t $OUTconc -u fasta -y raxml -c $CORES
echo concatenated finished

AMAS.py convert -i $OUTconc -f fasta -d $TYPE -u phylip -c $CORES
echo phylip format file obtained
