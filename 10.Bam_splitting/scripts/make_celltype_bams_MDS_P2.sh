#!/bin/bash

#SBATCH --job-name=MDS_P2_celltype_bams
#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=$1
metadata=$2
bam=$3

## step 1: generate list of barcodes for each celltype for MUT and WT 

cd $workdir 
mkdir barcode_lists
cd ./barcode_lists

#List of Cell Types in the file
sorted_unique_ids=($(cut -d" " -f7 $metadata | sort -u | tr '\n' ' '))


#cut -d " " -f 1,5,7: 1- barcode, 5- genotype, 7- celltype
for type in "${sorted_unique_ids[@]}"
do 
cut -d " " -f 1,5,7 $metadata | grep "WT $type$" | cut -d " " -f 1 | sed 's/\_2$//' > "$type"_wt_barcodes.txt 
cut -d " " -f 1,5,7 $metadata | grep "MUT $type$" | cut -d " " -f 1 | sed 's/\_2$//' > "$type"_mut_barcodes.txt 

mkdir "$type"
cd ./"$type"
mkdir split_bams
cd ./split_bams
mkdir wt_bams
mkdir mut_bams

## step 2: generate split bams for each cell type
cat "$workdir"/barcode_lists/"$type"_wt_barcodes.txt | while read barcode;
do
samtools view -h $bam | awk -v tag="BC:Z:$barcode" '($0 ~ /^@/ || index($0,tag)>0)' > wt_bams/"$barcode".sam
samtools view -S -b wt_bams/"$barcode".sam > wt_bams/"$barcode".bam
done

cat "$workdir"/barcode_lists/"$type"_mut_barcodes.txt | while read barcode;
do
samtools view -h $bam | awk -v tag="BC:Z:$barcode" '($0 ~ /^@/ || index($0,tag)>0)' > mut_bams/"$barcode".sam
samtools view -S -b mut_bams/"$barcode".sam > mut_bams/"$barcode".bam
done

cd ..
mkdir barcode_path_files
cd ./barcode_path_files

## step 3: merge split bams into one bam file 
awk 'NF{print "'$workdir'/barcode_lists/'$type'/split_bams/wt_bams/"$0".bam"}' "$workdir"/barcode_lists/"$type"_wt_barcodes.txt > "$type"_wt_bam_paths.csv

awk 'NF{print "'$workdir'/barcode_lists/'$type'/split_bams/mut_bams/"$0".bam"}' "$workdir"/barcode_lists/"$type"_mut_barcodes.txt > "$type"_mut_bam_paths.csv

cd ..
mkdir merge_bams
cd ./merge_bams

samtools merge "$type"_wt_merge.bam -b "$workdir"/barcode_lists/"$type"/barcode_path_files/"$type"_wt_bam_paths.csv
samtools merge "$type"_mut_merge.bam -b "$workdir"/barcode_lists/"$type"/barcode_path_files/"$type"_mut_bam_paths.csv

samtools sort "$type"_wt_merge.bam
samtools index "$type"_wt_merge.bam
samtools sort "$type"_mut_merge.bam
samtools index "$type"_mut_merge.bam

cd $workdir/barcode_lists
done

