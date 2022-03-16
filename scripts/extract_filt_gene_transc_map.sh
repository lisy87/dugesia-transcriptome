## script name: extract_filt_gene_transc_map.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N transc_map_filt             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"


##RUN FROM THE GENE_TRANSC_MAP FOLDER

PATH_names=$1
PATH_out=$2

for FILE in *.gene_trans_map
do
  SAMPLE=${FILE%%.gene_trans_map}""
  
  grep -wFf $PATH_names/$SAMPLE.filtered.txt $FILE > $PATH_out/$SAMPLE.gene_trans_map.filtered
    
done

