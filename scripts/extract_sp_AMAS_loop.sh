#script name: extract_sp_AMAS_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N ext_sp                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## activate conda

LIST=$1
TYPE=$2 #(dna,aa)
CORES=$3

SAMPLES=$(cat $LIST)

for file in *.fasta
do
	echo "*********************" extracting samples from $file "*************************"

	AMAS.py remove -x $SAMPLES -d $TYPE -f fasta -i $file -u fasta -c $CORES

	echo "********************" samples extracted from $file "***************************"

done
