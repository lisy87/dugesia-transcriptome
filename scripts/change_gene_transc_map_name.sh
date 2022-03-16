## script name: change_gene_transc_map_name.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N Gene_Transc_Map            #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

mkdir Gene_Transc_Map

for FOLDER in *.trinity_out
  do
  cd $FOLDER
  SAMPLE=${FOLDER%%.trinity_out}""
  mv Trinity.fasta.gene_trans_map $SAMPLE.gene_trans_map
  cp $SAMPLE.gene_trans_map ../Gene_Transc_Map
  cd ../

 done

