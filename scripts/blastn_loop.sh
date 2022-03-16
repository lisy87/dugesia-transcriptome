## script name: blastn_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N blastn             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate base

THREADS=$1


for FOLDER in *.BlobT

do

	SAMPLE=${FOLDER%%.BlobT}""

	cd $FOLDER

	echo "***********************" Blastn $SAMPLE "****************************"
	
  blastn -db PATH_to_nucleotide_database/nt -query *.fasta -num_threads $THREADS -max_target_seqs 10 -max_hsps 1 -outfmt '6 qseqid staxids bitscore std' -evalue 1e-25 > $SAMPLE.blastn.outfmt6

 	echo "*************************" Blastn $SAMPLE finish "**************************" 
  cd ../
done

