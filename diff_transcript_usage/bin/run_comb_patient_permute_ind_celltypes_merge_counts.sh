#!/bin/bash

#SBATCH --mem=64g
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/stdout_%j.log 

#########################################
### This script submits the script that permutes the combined patient data by merging all counts before calculating log odds ratio 
### permutations are done within each patient for each cell type matrix that is submitted to give one output 

#module load R/3.6.0

splitfiles=$1
genotype=$2
nperm=$3
patient=$4
outdir=$5
outfile=$6
runfiles=$7

echo $splitfiles
echo $genotype


Rscript "$runfiles"/combined_patient_permute_merge_reads_ind_celltypes_1WT.R $splitfiles $genotype $nperm $patient $outdir $outfile
