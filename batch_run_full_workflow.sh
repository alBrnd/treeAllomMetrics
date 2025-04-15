#!/bin/bash


echo "running part 1 of R script"
Rscript code/metrics_full_workflow_batch1.R


echo "running raycloudtools QSM"
bash code/batch_run_raycloudQSMs.sh

echo "running part 2 of R script"
Rscript code/metrics_full_workflow_batch2.R


echo "All scripts have been executed!"
