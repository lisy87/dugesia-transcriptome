## script name: assembly_stats.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N assemb_stats             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

#conda activate

printf filename'\t'total_length'\t'number'\t'mean_length'\t'longest'\t'shortest'\t'N_count'\t'Gaps'\t'N50'\t'N50n'\t'N70'\t'N70n'\t'N90'\t'N90n'\n' > assemblies_stats.tsv

for FILE in *.fasta

  do

  echo "*********************** ASSEMBLY-STATS" $FILE "*******************************************"

  assembly-stats -t -u $FILE >> assemblies_stats.tsv

  echo "*********************** ASSEMBLY-STATS" $FILE "FINISHED ********************************"

  done
