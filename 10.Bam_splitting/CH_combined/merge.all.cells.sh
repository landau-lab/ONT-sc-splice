#!/bin/bash

#SBATCH --job-name=CH_combined.merge.all
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

##Merge All WT bams
samtools merge CH_combined_all_merge.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/CH259_consensus.sorted.tags.GE.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/CH305_consensus.sorted.tags.GE.bam
samtools sort CH_combined_all_merge.bam
samtools index CH_combined_all_merge.bam


