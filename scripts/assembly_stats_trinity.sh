## script name: assembly_stats_trinity.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N assemb_stats_tr             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

mkdir Assembly_stats_Trinity

for FILE in *.fasta 
  do
    OUT=${FILE%.fasta}""
    echo "**************" Evaluating the assembly quality $FILE "***********************" >> Assembly_stats_Trinity/$OUT.assemb_stats.tsv

  PATH_to/trinityrnaseq-v2.9.1/util/TrinityStats.pl $FILE >> Assembly_stats_Trinity/$OUT.assemb_stats.tsv

done
