##script_name: FastQC.sh

##Put all the fasta files in the same folder and run from this directory

#!/bin/bash
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N fastqc             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate

mkdir FASTQC_results

for FILE in  *.gz
do
  echo "********************" $FILE "***************************"
  echo "*********************************************************"
  fastqc $FILE -o FASTQC_results/
done

# Run multiqc to obtain the single report of all samples
cd FASTQC_results/
multiqc .
