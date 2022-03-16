
#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N trimp_stats             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"


printf file'\t'seq_counts'\n' > trimp_stats.tsv

for f1 in *_1_trimP.fq.gz
do
  f2=${f1%%_1_trimP.fq.gz}""
  count=$(zcat $f1 | grep "@" | wc -l)
  printf $f2'\t'$count'\n' >> trimp_stats.tsv
done
