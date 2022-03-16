#script_name: trimAl.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N trimAl                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate

for FILE in *.fa
do
OUT=${FILE%%alig_oneline.trimmed.fa}"alig_oneline_trimmed_trimal.fasta"

trimal -in $FILE -out $OUT -automated1 -htmlout $OUT.html -keepheader

echo "*********************************************************************************"
echo "*************************" trimAl for $FILE done "*******************************"

done

