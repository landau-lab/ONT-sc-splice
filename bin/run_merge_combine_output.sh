#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=32g 

module load R/3.6.0

run_files=$1
outputsdir=$2
metadata=$3
finalout=$4

Rscript "$run_files"/merge_final_output_comb_patient_merge_counts.R $outputsdir $metadata $finalout


