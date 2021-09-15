#!/bin/bash

#SBATCH --job-name=CH259.merge.1UMI.all.mut.all.wt.bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=/gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH259

##Merge All WT bams
samtools merge CH259_all_wt_merge.1UMI.bam -b "$workdir"/all.wt.bams.paths.1UMI.csv
samtools sort CH259_all_wt_merge.1UMI.bam
samtools index CH259_all_wt_merge.1UMI.bam

##Merge All MUT bams
samtools merge CH259_all_mut_merge.1UMI.bam -b "$workdir"/all.mut.bams.paths.1UMI.csv
samtools sort CH259_all_mut_merge.1UMI.bam
samtools index CH259_all_mut_merge.1UMI.bam

