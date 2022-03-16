## script name: extract_filt_transc_BlobT.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N extract_transc             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate 

printf SAMPLE'\t'Platyhelminthes'\t'No-hit'\n' > blobtools_stats.tsv

for FOLDER in *.BlobT

do

	SAMPLE=${FOLDER%%.BlobT}""

	cd $FOLDER
  
  
  echo "*******************" stats Blobtools $SAMPLE "*************************"
	
	Platy=$(cat *.blobDB.table.txt | grep "Platyhelminthes" | wc -l)
	
	nohit=$(cat *.blobDB.table.txt | grep "no-hit" | wc -l) 
	
	printf $SAMPLE'\t'$Platy'\t'$nohit'\n' >> ../blobtools_stats.tsv
	
  
	echo "***********" Extract filtered transcripts $SAMPLE "*************"
	
	echo extracting Platyhelminthes
	cat *.blobDB.table.txt | grep "Platyhelminthes" | cut -f 1 > $SAMPLE.filtered.txt
	
	echo extracting no-hit
	cat *.blobDB.table.txt | grep "no-hit" | cut -f 1 >> $SAMPLE.filtered.txt
	
  
	echo changing transcripts name on the assembly
	cat *.fasta | cut -d " " -f 1 > $SAMPLE.newnames.fasta
	
	echo selecting filtered transcripts
	seqtk subseq $SAMPLE.newnames.fasta  $SAMPLE.filtered.txt > $SAMPLE.filt.transc.fasta
	
  echo "*****************" Filtered Transcripts Extracted $SAMPLE "*******************"
  cd ../
  
done

