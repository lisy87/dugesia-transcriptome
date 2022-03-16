#script name: view_folder_results.sh

#!/bin/sh

FOLDER=$1

cd $FOLDER
 echo files in $FOLDER
 ls -lrth | grep ".trinity_out"

 echo files in trinity_out folders
  for out in *.trinity_out
  do
      cd $out
      echo files in $out
      ls -lrth | grep "Trinity"
      cd ../
 done

cd ../
