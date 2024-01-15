#!/bin/bash

#SBATCH --mem=64g
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/stdout_%j.log 

#module load R/3.6.0

splitdir=$1
genotype=$2
nperm=$3
pattern=$4
outputdir=$5
outputfile=$6
runfiles=$7

echo $outputdir

Rscript "$runfiles"/split_JuncPermute_LogOR_perm_within_celltype_5p_3p.r \
  --split $splitdir \
  --genotype_file $genotype \
  --num_perm $nperm \
  --pattern $pattern \
  --output_dir $outputdir \
  --output_file $outputfile
