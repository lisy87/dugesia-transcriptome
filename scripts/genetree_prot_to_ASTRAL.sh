##script_name: genetree_prot_to_ASTRAL.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N iqtree_gene             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate iqtree

THREADS=$1

for FILE in *.fasta
do

TREE=${FILE%_prot.fasta}"_st.tree" #be sure to have the same format in the names of the files

iqtree -s $FILE -m LG+C20+F+G -ft $TREE -bb 10000 -nt $THREADS

echo "******************************" iqtree $FILE finished "****************************************"

done
