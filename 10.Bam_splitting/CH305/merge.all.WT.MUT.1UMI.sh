#!/bin/bash

#SBATCH --job-name=CH305.merge.all.mut.all.wt_bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=/gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH305

##Merge All WT bams
samtools merge CH305_all_wt_merge.bam -b "$workdir"/all.wt.bams.paths.csv
samtools sort CH305_all_wt_merge.bam
samtools index CH305_all_wt_merge.bam

##Merge All MUT bams
samtools merge CH305_all_mut_merge.bam -b "$workdir"/all.mut.bams.paths.csv
samtools sort CH305_all_mut_merge.bam
samtools index CH305_all_mut_merge.bam

