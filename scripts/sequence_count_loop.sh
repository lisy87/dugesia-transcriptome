## script name: sequence_count_loop.sh

#!/bin/bash

type=$1

for file in *.$type

do

	if [[ $type == *fastq* ]]
	then
		seqs=$(cat $file | grep "@" | wc -l)
		echo $seqs sequences in $file

	else
		seqs=$(cat $file | grep ">" | wc -l)
		echo $seqs sequences in $file

	fi

done

# Usage:  bash sequence_count_loop.sh fastq
