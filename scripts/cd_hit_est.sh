## script name: cd_hit_est.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N cd_hit             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate

## From cdhit help
## -c	sequence identity threshold, default 0.9
 	    #this is the default cd-hit's "global sequence identity" calculated as:
 	    #number of identical amino acids in alignment divided by the full length of the shorter sequence
## -G	use global sequence identity, default 1
	    # if set to 0, then use local sequence identity, calculated as :
 	    # number of identical amino acids in alignment divided by the length of the alignment
    	# NOTE!!! don't use -G 0 unless you use alignment coverage controls see options -aL, -AL, -aS, -AS
## -M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
## -T	number of threads, default 1; with 0, all CPUs will be used
## -n	word_length, default 10, see user's guide for choosing it
## -d	length of description in .clstr file, default 20
 	    #if set to 0, it takes the fasta defline and stops at first space
## -aL	alignment coverage for the longer sequence, default 0.0
 	    #if set to 0.9, the alignment must covers 90% of the sequence
## -aS	alignment coverage for the shorter sequence, default 0.0
 	    #if set to 0.9, the alignment must covers 90% of the sequence
## -g	1 or 0, default 0
 	    #by cd-hit's default algorithm, a sequence is clustered to the first cluster that meet the threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster that meet the threshold (accurate but slow mode)	but either 1 or 0 won't change the representatives of final clusters

echo creating CDHIT directory
mkdir CDHIT

THREADS=$1

for FILE in *.fasta
do
  OUT=${FILE%.Tr.fasta}"_Tr_cdhit.fasta"
  echo "*****************" running cdhit-est $FILE "*********************************************"
  cd-hit-est -i $FILE -o ./CDHIT/$OUT -c 0.99 -G 1 -M 120000 -T $THREADS -n 10 -d 0 -g 1  
  echo "******************************" cdhit-est $FILE finished "*********************************"
  echo "*******************************************************************************************"
done
