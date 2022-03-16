## script name: look_Dunuc12_in_OG.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N dunuc_OG    #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

SAMPLE=$1

TRANSC=$2


echo "**************" looking Dunuc12 for $SAMPLE "*******************"

	printf OG'\t'count'\n'> Dunuc12_all_OG.$SAMPLE.tsv

		for file in *.fa
		do

			OG=${file%%.fa}""

			seqcounts=$(cat $file | grep $TRANSC | wc -l)

			printf $OG'\t'$seqcounts'\n' >> Dunuc12_all_OG.$SAMPLE.tsv

		done

	awk '$2 != 0' Dunuc12_all_OG.$SAMPLE.tsv > Dunuc12_select_OG.$SAMPLE.tsv

	selectOG=$(cat Dunuc12_select_OG.$SAMPLE.tsv | wc -l)

echo "**************" Dunuc12 for $SAMPLE finish "*******************"
echo $selectOG OG with Dunuc12 for $SAMPLE
echo "***************************************************************"
