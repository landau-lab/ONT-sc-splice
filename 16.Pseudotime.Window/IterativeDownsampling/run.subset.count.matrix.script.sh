#!/bin/bash

#SBATCH --job-name=subset.count.matrix
#SBATCH --partition=bigmem,pe2
#SBATCH --mem=100g
#SBATCH --ntasks=10
#SBATCH --output=log.files/slurm-%j.out

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots

Rscript subset.count.matrix.R

