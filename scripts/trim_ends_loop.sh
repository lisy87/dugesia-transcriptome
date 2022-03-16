##script_name: trim_ends_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N trim_nuc                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

TRIMM_SCRIPT=$1 # PATH to script (for protein or nucleotide)

for file in *.fasta
do
	bash $TRIMM_SCRIPT $file
done
