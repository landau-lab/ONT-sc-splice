#!/bin/bash

#SBATCH --job-name=splice_annotator
#SBATCH --partition=bigmem
#SBATCH --mail-type=NONE
#SBATCH --cpus-per-task=12
#SBATCH --mem=200gb
#SBATCH --output=p6_10x_splice_annotator_stdout_%j.log

module load python/3.5.1

#python splice_annotator_got.py <path to datafiles> <input datafile> <junction tag .bed file> <output main file nam

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p1_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp1_10x_introns_primary_annotations_strand_filtered.txt 

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p2_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp2_10x_introns_primary_annotations_strand_filtered.txt 

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p3_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp3_10x_introns_primary_annotations_strand_filtered.txt

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data ch259_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed ch259_10x_introns_primary_annotations_strand_filtered.txt

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data ch305_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed ch305_10x_introns_primary_annotations_strand_filtered.txt

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p4_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp4_10x_introns_primary_annotations_strand_filtered.txt

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p5-1_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp5-1_10x_introns_primary_annotations_strand_filtered.txt 

#python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p5-2_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp5-2_10x_introns_primary_annotations_strand_filtered.txt  

python new_annotator_10x.py /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/3.10x_juncs/primary_annotations/data p6_10x_gtf_basic_strand_filtered.introns.info.txt leafviz_all_introns.bed mdsp6_10x_introns_primary_annotations_strand_filtered.txt 
