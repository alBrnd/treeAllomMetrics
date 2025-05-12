#!/bin/bash

data_folder="./data/rct_qsm/trees"

for data_file in "$data_folder"/*_raycloud.ply; do
  # Extract the basename without '_raycloud.ply'
  base_name=$(basename "$data_file" | sed 's/_raycloud\.ply//')

  # Find the corresponding mesh file
  mesh_file="$data_folder/${base_name}_raycloud_mesh.ply"
  
  # Check if the segmented file already exists
  segmented_file="$data_folder/${base_name}_raycloud_segmented.ply"
  if [ ! -f "$segmented_file" ]; then
    rayextract trees "$data_file" "$mesh_file" --max_diameter 2 --distance_limit 9 --height_min 5
  else
    echo "File $segmented_file already exists. Skipping rayextract."
  fi
done
