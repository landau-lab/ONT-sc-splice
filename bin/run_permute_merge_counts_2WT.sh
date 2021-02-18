#!/bin/bash

#SBATCH --mem=64g
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/stdout_%j.log 



module load R/3.6.0

splitfiles=$1
genotype=$2
nperm=$3
outdir=$4
outfile=$5
celltype=$6
runfiles=$7

echo $splitfiles
echo $genotype

Rscript "$runfiles"/combined_patient_permute_merge_reads_2WT_UMI.R $splitfiles $genotype $nperm $outdir $outfile $celltype
