#!/bin/bash

#SBATCH --job-name=quantile.counts
#SBATCH --partition=bigmem,pe2
#SBATCH --mem=300g
#SBATCH --ntasks=10
#SBATCH --output=log.files/slurm-%j.out

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots

#Rscript counts.per.pseudotime.quantile.MDS.R

#Rscript counts.per.pseudotime.quantile.CH259.R
#Rscript counts.per.pseudotime.quantile.CH305.R 

#Rscript counts.per.junc.MUTcells.only.MDS.962.R

#Rscript counts.per.junc.MUTcells.only.MDS.605.R

#Rscript counts.per.pseudotime.quantile.MDSP3.R

#Rscript counts.per.junc.MUTcells.only.CH305.R

#Rscript counts.per.junc.MUTcells.only.CH259.R

Rscript counts.per.junc.MUTcells.only.MDSP3.R
