#!/bin/bash

#SBATCH --job-name=celltype_bams
#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=$1
metadata=$2

cd $workdir/barcode_lists

sorted_unique_ids=($(cut -d" " -f8 $metadata | sort -u | tr '\n' ' '))

for type in "${sorted_unique_ids[@]}"
do
cut -d " " -f 1,6,8 $metadata | grep "WT $type$" | cut -d " " -f 1 | sed 's/\_1$//' > "$type"_wt_barcodes_2UMI.txt
cut -d " " -f 1,6,8 $metadata | grep "MUT $type$" | cut -d " " -f 1 | sed 's/\_1$//' > "$type"_mut_barcodes_2UMI.txt

cd ./"$type"/barcode_path_files

## step 2: Grab all the mut and wt cell barcodes that pass the threshold per cell type
awk 'NF{print "'$workdir'/barcode_lists/'$type'/split_bams/wt_bams/"$0".bam"}' "$workdir"/barcode_lists/"$type"_wt_barcodes_2UMI.txt > "$type"_wt_bam_paths_2UMI.csv

awk 'NF{print "'$workdir'/barcode_lists/'$type'/split_bams/mut_bams/"$0".bam"}' "$workdir"/barcode_lists/"$type"_mut_barcodes_2UMI.txt > "$type"_mut_bam_paths_2UMI.csv

cd ..

## step 3: merge split bams into one bam file 

cd ./merge_bams

samtools merge "$type"_wt_merge_2UMI.bam -b "$workdir"/barcode_lists/"$type"/barcode_path_files/"$type"_wt_bam_paths_2UMI.csv
samtools merge "$type"_mut_merge_2UMI.bam -b "$workdir"/barcode_lists/"$type"/barcode_path_files/"$type"_mut_bam_paths_2UMI.csv

samtools sort "$type"_wt_merge_2UMI.bam
samtools index "$type"_wt_merge_2UMI.bam
samtools sort "$type"_mut_merge_2UMI.bam
samtools index "$type"_mut_merge_2UMI.bam

cd $workdir/barcode_lists
done
