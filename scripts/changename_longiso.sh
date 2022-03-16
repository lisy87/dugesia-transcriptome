#script_name: changename_longiso.sh

#!/bin/bash

for FILE in *.fasta
do
  SAMPLE=${FILE%_longiso.fasta}""
  
  cat $FILE | cut -d " " -f 1 | sed 's/.p1//g' | sed "s/TRINITY/$SAMPLE/g" > $SAMPLE.longiso.newname.fasta 

done
