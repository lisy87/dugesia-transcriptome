##script_name: orthofinder_protein.sh

## Run into a folder with all files

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N ortfind_pt                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

#conda activate orthofinder3.6

THREADS=$1

PATH_to/OrthoFinder/orthofinder -f ./ -t $THREADS -M msa -S diamond_ultra_sens

echo "********************************************************************************************"
echo "*************************" ORTHOFINDER FINISHED "*******************************************"
