# script_name: change_name_SCnuc.sh 

#!/bin/sh

for FILE in *.fasta
do
        while read line
        do
        sed -i "s/$line//" $FILE
        done < species_listOK.txt
done

# remember to use the correct path to species_listOK.txt file

# For 1 file:
# while read line; do sed -i "s/$line//" OG0008560_nuc.fasta; done < species_listOK.txt

