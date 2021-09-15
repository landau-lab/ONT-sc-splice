#!/bin/bash
#SBATCH --job-name=samtools_split
#SBATCH --partition=pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=64gb
#SBATCH --output=logs/stdout_%j.log
#SBATCH --ntasks=1


source /etc/profile.d/modules.sh
module load java
module load samtools/1.10
module load bedtools


for i in $(cat ./4.Barcode_files/p1_barcodes_trimmed.tsv); do

#i=$1
#samtools view -d BC:Z:${i} -o /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/split_by_bc/"$i".bam ./1.Bam_files/1.ONT_sicelore_output/p1_consensus.tags.GE_no_secondary_sorted.bam

samtools view -h -F 256 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i".bam > /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i"_temp.bam

samtools merge -h /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_header.sam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i"_filt.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i"_temp.bam /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/empty.bam

rm /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i"_temp.bam

bedtools bamtofastq -i /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i"_filt.bam -fq /gpfs/commons/groups/landau_lab/SF3B1_splice_project/1.Bam_files/1.ONT_sicelore_output/p1_split_by_bc/"$i".fastq

done

