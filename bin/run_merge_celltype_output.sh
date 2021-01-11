#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=32g 

module load R/3.6.0

run_files=$1
outputsdir=$2
metadata=$3
genotype=$4
finalout=$5

Rscript "$run_files"/merge_ind_celltype_output.R $outputsdir $metadata $genotype $finalout


