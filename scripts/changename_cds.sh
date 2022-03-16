#script_name: changename_cds.sh

#!/bin/bash

for FILE in *.cds
do
  SAMPLE=${FILE%.filt.transc.fasta.transdecoder.cds}""
  
  cat $FILE | cut -d " " -f 1 | sed 's/.p1//g' | sed "s/TRINITY/$SAMPLE/g" > $SAMPLE.cds.newname.fasta 
done

