#!/bin/bash

#SBATCH --job-name=CH_combined.merge.all.mut.all.wt_bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

##Merge All WT bams
samtools merge CH_combined_all_wt_merge.2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH305/CH305_all_wt_merge.2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH259/CH259_all_wt_merge.2UMI.bam
samtools sort CH_combined_all_wt_merge.2UMI.bam
samtools index CH_combined_all_wt_merge.2UMI.bam

##Merge All MUT bams
samtools merge CH_combined_all_mut_merge.2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH305/CH305_all_mut_merge.2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH259/CH259_all_mut_merge.2UMI.bam
samtools sort CH_combined_all_mut_merge.2UMI.bam
samtools index CH_combined_all_mut_merge.2UMI.bam

