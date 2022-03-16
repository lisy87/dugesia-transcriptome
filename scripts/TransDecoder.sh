## script name: TransDecoder.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N transDeco             #name Job
#$ -m ea                 #add email
echo "************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "************************************************************"

#conda activate transdecoder3.6

THREADS=$1
PATH_Gen_transc_map=$2
 
for FILE in *.fasta
do

  SAMPLE=${FILE%%.filt.transc.fasta}""
  
  echo "******************************************************************"
  echo "*******************" transdecoder $FILE "*************************"

  TransDecoder.LongOrfs -t $FILE --gene_trans_map $PATH_Gen_transc_map/$SAMPLE.gene_trans_map.filtered

  TransDecoder.Predict -t $FILE

done


