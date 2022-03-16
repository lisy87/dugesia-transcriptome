##script_name: longiso_count_loop.sh

#!/bin/bash

type=$1
OUT=$2
for file in *.$type

do

	if [[ $type == *fastq* ]]
	then
		seqs=$(cat $file | grep "@" | wc -l)
		echo $seqs sequences in $file >> $OUT

	else
		seqs=$(cat $file | grep ">" | wc -l)
		echo $seqs sequences in $file >> $OUT

	fi

done

##Usage: In a folder with fasta (or fastq) files run: bash longiso_count_loop.sh fasta/fastq out.txt

