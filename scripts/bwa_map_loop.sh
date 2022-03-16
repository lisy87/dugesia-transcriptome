## script name: bwa_map_loop.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N bwa_map             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

#conda activate snake3.6

THREADS=$1
PATHtoreads=$2

for FOLDER in *.BlobT

do

	SAMPLE=${FOLDER%%.BlobT}""

	READS1=${FOLDER%%.BlobT}"_1_trimP.fq.gz"

	READS2=${FOLDER%%.BlobT}"_2_trimP.fq.gz"

	cd $FOLDER

	echo "***********************" Indexing $SAMPLE "****************************"

	bwa index *.fasta

	echo "*************************" Mapping $SAMPLE "**************************"

	bwa mem -t $THREADS -o $SAMPLE.sam *.fasta  $PATHtoreads/$READS1  $PATHtoreads/$READS2

	echo "************************" SAM to BAM "*********************************"

	samtools view -S $SAMPLE.sam -bo $SAMPLE.bam

	echo "***********************" sorting BAM "*********************************"

	samtools sort -o $SAMPLE.sorted.bam $SAMPLE.bam

	echo "*******************" indexing BAM "************************************"

	samtools index $SAMPLE.sorted.bam
 
	echo "*******************" deleting intermediate files "*********************"

	rm $SAMPLE.sam
	rm $SAMPLE.bam

	echo "************************" $SAMPLE Mapped " **************************"

	echo "************************" Stats $SAMPLE Mapping "********************"

	echo "************************" Stats $SAMPLE Mapping "********************" >> ../Bam_stats.txt
	samtools flagstat $SAMPLE.sorted.bam >> ../Bam_stats.txt
	echo "*********************************************************************" >> ../Bam_stats.txt

	cd ../
done

