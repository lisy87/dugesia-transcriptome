## script name: Trinity_assembly_remove_extrafiles.sh

##Put reads_1 and reads_2 in the same folder and run from this directory

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N Trinity             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate snake3.6

MEMORY=$1

CPU=$2

for READ1 in *_1_trimP.fq.gz

  do
    READ2=${READ1%%_1_trimP.fq.gz}"_2_trimP.fq.gz"

    NAME=${READ1%%_1_trimP.fq.gz}""

    echo "**************" Assembly $READ1 and $READ2 "*************"

Trinity --seqType fq --max_memory $MEMORY --left $READ1 --right $READ2 --CPU $CPU --output ./$NAME.trinity_out --full_cleanup

      echo "**************" $NAME Assembled "*************"
      echo "*******************************************************"
done

echo "*********************************************************************************************"
echo The assemblies are finished
