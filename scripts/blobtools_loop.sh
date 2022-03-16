## script name: blobtools_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N blobtools             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate blobtools3.6

for FOLDER in *.BlobT

do

	SAMPLE=${FOLDER%%.BlobT}""
  FASTAFILE=${FOLDER%%.BlobT}"_Tr_cdhit.fasta"
	cd $FOLDER

	echo "***********************" Blobtools $SAMPLE "****************************"
	
	echo Creating blobDB
	
  PATH_to_desired_directory/blobtools create -i $FASTAFILE -b $SAMPLE.sorted.bam -x bestsumorder -t $SAMPLE.blastn.outfmt6 -o ./$SAMPLE.blobt.out
  
  echo Creating view
  
  PATH_to_desired_directory/blobtools view -i *.blobDB.json  -x bestsumorder
  
  echo Creating plot
  
  PATH_to_desired_directory/blobtools plot -i *.blobDB.json -x bestsumorder
  

  echo "***********************" Blobtools $SAMPLE finish "****************************"
  cd ../
  
done

