## script name: Trinity_assembly_remove_extrafiles.sh

##Put reads_1 and reads_2 in the same folder and run from this directory

#!/bin/sh
#$ -cwd
#$ -j y
#$ -V                    #export environment var
#$ -N Trinity             #name Job
#$ -m ea                 #add email
echo "***************************************************************"
echo "*********" $HOSTNAME " ****** JOB_ID=" $JOB_ID "  *************"
echo "***************************************************************"

## conda activate snake3.6

MEMORY=$1

CPU=$2

for READ1 in *_1_trimP.fq.gz

  do
    READ2=${READ1%%_1_trimP.fq.gz}"_2_trimP.fq.gz"

    NAME=${READ1%%_1_trimP.fq.gz}""

    echo "**************" Assembly $READ1 and $READ2 "*************"

Trinity --seqType fq --max_memory $MEMORY --left $READ1 --right $READ2 --CPU $CPU --output ./$NAME.trinity_out

      echo "**************" $NAME Assembled "*************"
      echo "**************" removing extra files "******************"
        cd ./$NAME.trinity_out
        rm left.fa.ok
        rm right.fa.ok
        rm both.fa
        rm both.fa.ok
        rm both.fa.read_count
        rm jellyfish.kmers.25.asm.fa
        rm jellyfish.kmers.25.asm.fa.histo
        rm inchworm.kmer_count
        rm inchworm.DS.fa
        rm inchworm.DS.fa.finished
        rm scaffolding_entries.sam
        rm pipeliner.17805.cmds
        rm -r chrysalis
        rm read_partitions
        rm partitioned_reads.files.list
        rm partitioned_reads.files.list.ok
        rm recursive_trinity.cmds
        rm recursive_trinity.cmds.ok
        rm recursive_trinity.cmds.completed
        rm -r insilico_read_normalization
        cd ../

      echo "**************" extra files removed "******************"
      echo "*******************************************************"
done

echo "*********************************************************************************************"
echo The assemblies are finished
