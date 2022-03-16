#script_name: samenameTOconcat.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N samename                   #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

for file in *.fasta
do
  cat $file | sed 's/_DN/@/' | cut -d "@" -f 1 > $file.samenames
done

# Use filextention.sh to change the names of the outputs 

