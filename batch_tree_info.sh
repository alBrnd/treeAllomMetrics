#!/bin/bash

for file in ./rct_qsm/trees/*trees.txt; do
  if [ -f "$file" ]; then
    echo "Processing $file"
    # Add your commands here
    treeinfo "$file" --branch_data
 fi
done
