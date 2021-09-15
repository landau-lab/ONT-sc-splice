#!/bin/bash

#SBATCH --job-name=celltype_bams
#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2


module load samtools

workdir=$1
metadata=$2
gt_field=$3
ct_field=$4
bam=$5

## step 1: generate list of barcodes for each celltype for MUT and WT 

#cd $workdir
#mkdir barcode_lists
#cd ./barcode_lists

#List of Cell Types in the file
echo $ct_field

sorted_unique_ids=($(cut -d" " -f "$ct_field" $metadata | sort -u | tr '\n' ' '))

for type in "${sorted_unique_ids[@]}"
do
   echo $type

done


#sorted_unique_ids=($(cut -d" " -f10 mds_p1.genotype.info.with.cell.type.txt | sort -u | tr '\n' ' '))
#echo $sorted_unique_ids

# Read the array values with space
#for val in "${sorted_unique_ids[@]}"; do
#  echo $val
#done
