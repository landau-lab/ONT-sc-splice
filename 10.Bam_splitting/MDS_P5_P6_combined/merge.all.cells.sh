#!/bin/bash

#SBATCH --job-name=MDS.untreated_combined.merge.all
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

##Merge All WT bams
samtools merge MDS_untreated_combined_all_merge.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P5_1_consensus.sorted.tags.GE.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P5_2_consensus.sorted.tags.GE.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P6_consensus.sorted.tags.GE.bam
samtools sort MDS_untreated_combined_all_merge.bam
samtools index MDS_untreated_combined_all_merge.bam


