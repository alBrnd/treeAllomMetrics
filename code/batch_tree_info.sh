#!/bin/bash



for file in ./data/rct_qsm/trees/*trees.txt; do
  # Extract the basename without '.txt
  base_name=$(basename "$file" | sed 's/\.txt//')

  # Check if the info file already exists
  info_file="./data/rct_qsm/trees/${base_name}_info.txt"	

  if [ ! -f "$info_file" ]; then
    echo "Processing $file"
    treeinfo "$file" --branch_data
  else
    echo "File $info_file already exists. Skipping info."	
  fi
done

