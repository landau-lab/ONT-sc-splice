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

Rscript "$runfiles"/split_JuncPermute_LogOR_combined_patient_within_celltype_MDS_treated_remove_celltype.R $splitfiles $genotype $nperm $outdir $outfile $celltype
