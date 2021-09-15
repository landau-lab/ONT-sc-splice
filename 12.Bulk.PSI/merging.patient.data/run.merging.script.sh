#!/bin/bash
#SBATCH --job-name=merging.count.matricies
#SBATCH --partition=bigmem,pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=400g

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots

#Rscript merging.count.matricies.R count.matricies.paths.CH.txt CH259.CH305

Rscript merging.count.matricies.R count.matricies.paths.MDS.txt MDSP5.MDSP6
