
##script_name: ASTRALpro_default.sh

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N ASTRALdef             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

# conda activate astral

INPUT_TREES=$1 #all gene trees in one file all_gene_trees.tree
THREADS=$2
OUT=$3
LOG=$4

java -D"java.library.path=PATH/A-o-master/ASTRAL-MP/lib" -jar PATH/A-o-master/ASTRAL-MP/astral.x.x.x.jar -i $INPUT_TREES -T $THREADS -o $OUT 2>$LOG

# NOTE: Use your own PATH to A-o-master/ASTRAL-MP/ and substitute astral.x.x.x.jar by the installed version
