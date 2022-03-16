## script name: sequences_count.sh

#!/bin/bash

file=$1
type=$2

if [[ $type == *fastq* ]]
then
	seqs=$(cat $file | grep "@" | wc -l)
	echo $seqs sequences in $file

else
	seqs=$(cat $file | grep ">" | wc -l)
	echo $seqs sequences in $file

fi

# Usage: bash ./sequence_count.sh example.fasta fasta
