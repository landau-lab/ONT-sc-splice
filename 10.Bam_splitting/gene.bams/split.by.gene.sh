#!/usr/bin/env bash

##Format for command line run: bash split.by.gene.sh path/to/bamfile gene_name gene_coordinates patientID
##Example: bash split.by.gene.sh /gpfs/commons/home/pchamely/robbie/ont_bams PIK3CA chr3:179146357-179242093 CH

module load samtools

bamfile=$1
gene_name=$2
chr_region=$3
patient_id=$4

samtools view -h "$bamfile" "$chr_region" > "$patient_id"_"$gene_name".bam
samtools sort "$patient_id"_"$gene_name".bam > "$patient_id"_"$gene_name"_sorted.bam
samtools index "$patient_id"_"$gene_name"_sorted.bam
