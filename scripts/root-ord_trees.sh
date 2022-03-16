## script_name: root-ord_trees.sh

#!/bin/bash

#conda activate

for file in *.tree
do
nw_reroot $file -l "outgroup" | nw_order - > $file.rooted.order
echo "***********" $file rooted "*************"
done
