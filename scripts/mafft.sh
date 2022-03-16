##script_name: mafft.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N mafft                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

##activate conda

THREADS=$1

mkdir mafft_alignments

for FILE in *.fasta
do
  GENE=${FILE%.fasta}""
  OUT=${FILE%.fasta}"_alig.fasta"
  mafft --auto --maxiterate 1000 --thread $THREADS $FILE > ./mafft_alignments/$OUT
  echo "*********************" $GENE aligned "************************"
done

