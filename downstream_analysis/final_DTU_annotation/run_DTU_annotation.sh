#!/bin/bash

#SBATCH --mem=32g
#SBATCH --partition=pe2
#SBATCH --job-name=submit_annotation

module load R/3.6.0

threedata=$1
fivedata=$2
output=$3

Rscript /gpfs/commons/groups/landau_lab/mariela/tools/ONT-sc-splice/bin/DTU_junction_annotation.R $threedata $fivedata $output

