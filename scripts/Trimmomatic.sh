##script_name: Trimmomatic.sh

##Put reads1 and reads2 in the same folder and run from said directory.
# The files have to have the following structure in the names:
# samplename_1.fastq.gz reads 1
# samplename_2.fastq.gz reads 2

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N Trimm             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate trimmo

THREADS=$1
ADAPTERS=$2 ## PATH to ADAPTERS

for f1 in *_1.fastq.gz
do
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
   trimmomatic PE -threads $THREADS -phred33 -trimlog $f1.log $f1 $f2  $f1.trimP.fq.gz $f1.trimUP.fq.gz $f2.trimP.fq.gz $f2.trimUP.fq.gz ILLUMINACLIP:$ADAPTERS:2:30:10:1:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
