#!/bin/bash

for file in ./rct_qsm/trees/*_raycloud.ply; do
  if [ -f "$file" ]; then
    echo "Processing $file"
    # Add your commands here
    rayextract terrain "$file"
  fi
done
