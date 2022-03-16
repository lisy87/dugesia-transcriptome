##script_name: Trimmomatic_uniq_index.sh

##Put reads1, reads2, and illuminaclimp files in the same folder and run from said directory. 
# The files have to have the following structure in the names:
# samplename_1.fastq.gz reads 1
# samplename_2.fastq.gz reads 2
# samplename_TruSeq_adap_ind.fa illuminaclimp file

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

for f1 in *_1.fastq.gz
do
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
    f3=${f1%%_1.fastq.gz}"_TruSeq_adap_ind.fa"
   trimmomatic PE -threads 8 -phred33 -trimlog $f1.log $f1 $f2  $f1.trimP.fq.gz $f1.trimUP.fq.gz $f2.trimP.fq.gz $f2.trimUP.fq.gz ILLUMINACLIP:$f3:2:30:10:1:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
