## script name: create_blobtools_folders.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N blobtFolder             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

mkdir BLOBTOOLS

for FILE in *_Tr_cdhit.fasta

do

  SAMPLE=${FILE%%_Tr_cdhit.fasta}""
  mkdir BLOBTOOLS/$SAMPLE.BlobT
  cp $FILE BLOBTOOLS/$SAMPLE.BlobT

done

