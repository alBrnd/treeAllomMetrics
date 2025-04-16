#!/bin/bash
rm data/rct_qsm/trees/*_raycloud*

#echo "Running txt2ply.py on all trees in directory..."
#python convert_to_ply.py trees trees

echo "Running rayimport..."
bash code/batch_rayimport.sh

echo "Running rayextract terrain..."
bash code/batch_rayextract_terrain.sh

echo "Running rayextract trees..."
bash code/batch_rayextract_trees.sh

echo "Running batch_tree_info.sh"
bash code/batch_tree_info.sh

echo "All scripts have been executed!"
