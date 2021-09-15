#!/bin/bash
#SBATCH --job-name=bulk.PSI
#SBATCH --partition=bigmem,pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=200g
#SBATCH --output=err.and.out/CH259.CH305.merged_stdout_%j.log

module load R/3.5.1

patientID=CH259.CH305.merged 

Rscript bulkPSI.indiv.patient.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/"$patientID".count.matrix.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/"$patientID".intron.metadata.txt output/"$patientID".all.intron.metadata.with.bulkPSI.ONT.tsv

