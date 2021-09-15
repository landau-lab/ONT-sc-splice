#!/bin/bash

#SBATCH --job-name=1UMI.MDS.untreated.combined.merge.all.mut.all.wt_bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

##Merge All WT bams
samtools merge MDS.untreated_combined_all_wt_merge.1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1/bulk.bams/MDS_P5_1_all_wt_merge_1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_2/bulk.bams/MDS_P5_2_all_wt_merge_1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P6/bulk.bams/MDSP6_all_wt_merge_1UMI.bam
samtools sort MDS.untreated_combined_all_wt_merge.1UMI.bam
samtools index MDS.untreated_combined_all_wt_merge.1UMI.bam

##Merge All MUT bams
samtools merge MDS.untreated_combined_all_mut_merge.1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1/bulk.bams/MDS_P5_1_all_mut_merge_1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_2/bulk.bams/MDS_P5_2_all_mut_merge_1UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P6/bulk.bams/MDSP6_all_mut_merge_1UMI.bam
samtools sort MDS.untreated_combined_all_mut_merge.1UMI.bam
samtools index MDS.untreated_combined_all_mut_merge.1UMI.bam
