##script_name: BUSCO_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N busco            #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate busco

CORES=$1

for file in *.fasta
do
  OUT=${file%%_Tr.fasta}"_busco"
  busco -i $file -l metazoa_odb10 -o $OUT -m transcriptome -c $CORES
done

