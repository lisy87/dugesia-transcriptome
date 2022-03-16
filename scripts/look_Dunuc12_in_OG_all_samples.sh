## script name: look_Dunuc12_in_OG_all_samples.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N dunuc_OG_all    #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

TRANSC=$1


echo "**************" looking Dunuc12 for ALL SAMPLEs "*******************"

	printf OG'\t'count'\n'> Dunuc12_all_OG.ALL_SAMPLES.tsv

		for file in *.fa
		do

			OG=${file%%.fa}""

			seqcounts=$(cat $file | egrep 'DetruBerg_2_DN1077_c0_g1 | DetruBerg_3_DN3226_c0_g1 | DetruBerg_5_DN6190_c0_g1 | DliguCat_3_DN1320_c0_g1 | DliguTri_C_DN596_c0_g1 | DliguTri_D_DN784_c0_g1 | DetruPie_2_DN92587_c0_g1 | DliguGarda_1_DN77494_c0_g1' | wc -l)

			printf $OG'\t'$seqcounts'\n' >> Dunuc12_all_OG.ALL_SAMPLES.tsv

		done

	awk '$2 != 0' Dunuc12_all_OG.ALL_SAMPLES.tsv > Dunuc12_select_OG.ALL_SAMPLES.tsv

	selectOG=$(cat Dunuc12_select_OG.ALL_SAMPLES.tsv | wc -l)

echo "**************" Dunuc12 for ALL_SAMPLES finish "*******************"
echo $selectOG OG with Dunuc12 for ALL_SAMPLES
echo "***************************************************************"
