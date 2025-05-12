#!/bin/bash


for file in ./data/rct_qsm/trees/*_raycloud.ply; do
  # Extract the basename without '.ply'
  base_name=$(basename "$file" | sed 's/\.ply//')

  # Check if the mesh file already exists
  mesh_file="./data/rct_qsm/trees/${base_name}_mesh.ply"	

  if [ ! -f "$mesh_file" ]; then
    echo "Processing $file"
    rayextract terrain "$file"
  else
    echo "File $mesh_file already exists. Skipping terrain."	
  fi
done

