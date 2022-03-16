#script name: create_folders.sh

#!/bin/sh

FOLDERS=$1 # mount of folders to create
FILES=$2 # mount of files to move to each folder

for i in `eval echo {1..$FOLDERS}`; do mkdir $i.FOLDER; done

for folder in *.FOLDER; do ls -1 *.fasta | head -n $FILES | xargs -I{} mv {} $folder; done

