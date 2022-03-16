##script_name: interleave_to_oneline_loop.sh

#!/bin/sh
for FILE in *.fasta
  do
  OUT=${FILE%%.fasta}"_oneline.fasta"
  perl PATH_to/interleave_to_oneline.pl $FILE > $OUT
done
