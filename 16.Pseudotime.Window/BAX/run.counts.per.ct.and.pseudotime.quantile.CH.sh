#!/bin/bash

#SBATCH --job-name=quantile.ct.counts
#SBATCH --partition=bigmem,pe2
#SBATCH --mem=100g
#SBATCH --ntasks=10
#SBATCH --output=slurm-%j.out

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots


Rscript counts.per.ct.and.pseudotime.quantile.CH259.305.R
