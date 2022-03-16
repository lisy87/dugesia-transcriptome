##script_name: Filextention_change.sh

## Put all the files in a folder and run: Filextention_change.sh file_extension_1  file_extension_2

#!/bin/bash

for filextention in *$1
  do
        mv $filextention ${filextention%$1}$2
done

# Usage: For change the file .fa to .fasta
# Filextention_change.sh .fa .fasta
