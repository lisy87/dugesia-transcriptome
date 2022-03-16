#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N Trinity_out             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

mkdir Trinity_fastas

for FOLDER in *.trinity_out
  do
  cd $FOLDER
  SAMPLE=${FOLDER%%.trinity_out}""
  mv Trinity.fasta $SAMPLE.Tr.fasta
  cp $SAMPLE.Tr.fasta ../Trinity_fastas
  cd ../

 done

