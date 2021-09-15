#!/bin/bash

#SBATCH --job-name=2UMI.MDS_untreated_combined_celltype_bams
#SBATCH --mem=200g
#SBATCH --partition=bigmem,pe2

module load samtools

celltypes=("HSPC" "IMP" "EP" "MEP" "NP" "MkP" "MonoDC" "Mat_Mono" "DC" "Bcells" "Tcells" "NA")


for type in ${celltypes[@]}
do

echo "$type"
#mkdir "$type"
cd ./"$type"

#Merge all the per patient WT bams 
samtools merge MDS_untreated_combined_"$type"_wt_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1/barcode_lists/"$type"/merge_bams/"$type"_wt_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_2/barcode_lists/"$type"/merge_bams/"$type"_wt_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P6/barcode_lists/"$type"/merge_bams/"$type"_wt_merge_2UMI.bam
  
samtools sort MDS_untreated_combined_"$type"_wt_merge_2UMI.bam
samtools index MDS_untreated_combined_"$type"_wt_merge_2UMI.bam

#Merge all the per patient MUT bams
samtools merge MDS_untreated_combined_"$type"_mut_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1/barcode_lists/"$type"/merge_bams/"$type"_mut_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_2/barcode_lists/"$type"/merge_bams/"$type"_mut_merge_2UMI.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P6/barcode_lists/"$type"/merge_bams/"$type"_mut_merge_2UMI.bam

samtools sort MDS_untreated_combined_"$type"_mut_merge_2UMI.bam
samtools index MDS_untreated_combined_"$type"_mut_merge_2UMI.bam

cd ..
done
