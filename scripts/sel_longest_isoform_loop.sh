##script_name: sel_longest_isoform_loop.sh

# Script to use the choose_longest_iso.py script of Tauana Cunha  http://dx.doi.org/10.1098/rspb.2018.2776 into a loop with multiple files of protein sequences obtained from Transdecoder analysis.  


#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N sel_isof             #name Job
#$ -m ea                 #add email
echo "************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "************************************************************"

##use python2.7 to run this script

for FILE in *.pep
  do
	OUT=${FILE%.pep}"_longiso.fasta"
	echo "**************************" selecting longest isoform in $FILE "****************************"
	python PATH_to_script/choose_longest_iso.py -l -i=$FILE -o=$OUT
	echo "********************************************************************************************"
	echo longest isoform selection finished
done

