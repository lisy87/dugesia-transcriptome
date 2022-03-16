## script name: pep_stats.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N pep_stats             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"


printf sample'\t'peptides'\t'complete'\t'missingstart'\t'missingstop'\t'incomplete'\n' > pepstats.tsv

for FILE in *.pep

  do

  sample=${FILE%.pep}""
  
  peptides=$(cat $FILE | grep ">" | wc -l)

  complete=$(cat $FILE | grep "complete" | wc -l)

  missingstart=$(cat $FILE | grep "5prime_partial" | wc -l)

  missingstop=$(cat $FILE | grep "3prime_partial" | wc -l)

  incomplete=$(cat $FILE | grep "internal" | wc -l)

  printf $sample'\t'$peptides'\t'$complete'\t'$missingstart'\t'$missingstop'\t'$incomplete'\n' >> pepstats.tsv

done

