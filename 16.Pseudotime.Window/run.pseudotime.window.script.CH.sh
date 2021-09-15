#!/bin/bash

#SBATCH --job-name=pseudotime.window
#SBATCH --partition=bigmem,pe2
#SBATCH --mem=300g
#SBATCH --ntasks=10
#SBATCH --output=log.files/CH_slurm-%j.out

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots

#Rscript PseudotimeWindow.CH.individual.R

#Rscript PseudotimeWindow.CH.merged.R  

Rscript PseudotimeWindow.CH.merged.onlyMUT.R

