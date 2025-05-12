#!/bin/bash

for file in ./data/rct_qsm/trees/*filtered.ply; do
  # Extract the basename without '.ply'
  base_name=$(basename "$file" | sed 's/\.ply//')

  # Check if the raycloud file already exists
  raycloud_file="./data/rct_qsm/trees/${base_name}_raycloud.ply"	

  if [ ! -f "$raycloud_file" ]; then
    echo "Processing $file"
    rayimport "$file" 0,0,-1 --max_intensity 0 --remove_start_pos
  else
    echo "File $raycloud_file already exists. Skipping rayimport."	
  fi
done
