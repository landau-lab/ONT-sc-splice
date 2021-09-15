#!/bin/bash

#SBATCH --job-name=MDS_P5_1.merge.all.mut.all.wt_bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=/gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1/bulk.bams


##Merge All WT bams
samtools merge MDS_P5_1_all_wt_merge_1UMI.bam -b "$workdir"/all.wt.bams.paths.1UMI.csv
samtools sort MDS_P5_1_all_wt_merge_1UMI.bam
samtools index MDS_P5_1_all_wt_merge_1UMI.bam

##Merge All MUT bams
samtools merge MDS_P5_1_all_mut_merge_1UMI.bam -b "$workdir"/all.mut.bams.paths.1UMI.csv
samtools sort MDS_P5_1_all_mut_merge_1UMI.bam
samtools index MDS_P5_1_all_mut_merge_1UMI.bam

