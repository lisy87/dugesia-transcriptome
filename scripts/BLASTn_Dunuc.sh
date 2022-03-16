# script name: BLASTn_Dunuc.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N blastn_dunuc             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

#conda activate base

THREADS=$1

TRANSC=$2

SAMPLE=${TRANSC%%.filt.transc.fasta}""

echo "***********************" Creating BLAST database of $SAMPLE "****************************"

makeblastdb -in $TRANSC -dbtype nucl -parse_seqids

echo "***********************" BLAST database of $SAMPLE created "****************************"
echo "****************************************************************************************"


for file in *_CS.fasta
do
        CLONE=${file%%_CS.fasta}""

echo "********************" Blasting $CLONE "****************************************"

blastn -db $TRANSC -query $file -num_threads $THREADS -max_target_seqs 10 -max_hsps 1 -outfmt '6 qseqid bitscore std ' > $CLONE.blastn.outfmt6

echo "********************" Blasting $CLONE finished "****************************************"
echo "****************************************************************************************"

done
