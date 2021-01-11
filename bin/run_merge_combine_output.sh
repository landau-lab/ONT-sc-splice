#!/bin/bash 

#SBATCH --partition=bigmem,pe2
#SBATCH --mem=100g
#SBATCH --cpus-per-task=4

module load R/3.6.0

runfiles=$1
resultsout=$2
comb_metadata=$3
patients=$4
output=$5

Rscript "$runfiles"/merge_final_output_comb_patient.R $resultsout $comb_metadata $patients $output
