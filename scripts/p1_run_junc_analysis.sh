#!/bin/bash

#SBATCH --job-name=P1_junc_analysis
#SBATCH --cpus-per-task=12
#SBATCH --mem=200gb
#SBATCH --partition=bigmem

source /etc/profile.d/modules.sh
module load R
module load anaconda3
module loady python
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/ahawkins/.conda/envs/bio

echo "Completing strand adjustment"
#arg1 is leafcutter counts matrix 
#arg2 is leafcutter metadata converted to csv
#arg 3 is directory for outputs - strand adjusted counts and strand adjusted matrix

		#-------- Run strand adjustment of matrix -----------#
#Rscript ./scripts/strand_adjustment.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/1.Original_Leafcutter_Output/p1_align_smartseq_nofilter_perind_numers.counts.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/1.Original_Leafcutter_Output/p1_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/2.Strand_Adjusted_Matrix/


echo "Strand adjustment done" 

echo "Filtering counts matrix" 
#arg1 is strand adjusted counts matrix
#arg2 is strand adjusted metadata 
#arg3 is output directory - counts filtered matrix and counts filtered metadata 
#arg4 is minimum number of reads to filter on 
		#--------- Filter counts matrix for total reads covering junction > 10 -------# 
Rscript ./scripts/filter_counts_matrix.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/2.Strand_Adjusted_Matrix/strand_adjusted_counts.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/2.Strand_Adjusted_Matrix/strand_adjusted_metadata.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/ 10


echo "Counts matrix filtered" 

                #----- copy input files for annotation to folder -----# 
cp /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/metadata.filtered.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/4.Leafcutter_Annotation/metadata.filtered.txt 

echo "Annotating metadata" 
		#--------- Adding annotations to metadata ------------#
#arg1 is directory to input files
#arg2 is filtered metadata matrix
#arg3 is gtf 
#arg4 is name of output of annotations
#arg5 is name of output of exon information 

python ./scripts/splice_annotator_got_v4.py /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/4.Leafcutter_Annotation/ metadata.filtered.txt gencode.v31.basic.annotation.csv metadata_annotations.csv exons.csv
echo "Annotations done"

echo "Running junction permutations" 
#arg1 is filtered counts matrix
#arg2 is genotyping information 
#arg3 is nperm (need to actually change in R script though) 
#arg4 is output directory - Rdata with pvalues and observed difference between MUT - WT normalized counts 

		#----------Run Mut vs. WT Junction Permutation Test-----------#
#100k permutations
#Rscript ./scripts/runJuncPermutation_mut_wt_v2.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/counts.filtered.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt 100000 /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/5.Junction_permutation/100k_perm_v3 


#Rscript ./scripts/runJunctionPermutation_mut_wt.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/P1_counts.filtered.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt 100000 /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/5.Junction_permutation/100k_perm/

echo "Calculating log OR"
#arg1 is filtered counts matrix 
#arg2 is genotyping info
#arg3 is filtered metadata
#arg4 is output directory - Rdata with pvalues and observed log odds ratio for each junction in compared to all other junctions in that cluster 
 
	#---------- Run Mut vs. WT Log OR for each junction with cluster > 1 -------#
#Rscript ./scripts/JuncPermute_LogOR.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/counts.filtered.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/4.Leafcutter_Annotation/P1_metadata_annotations.csv /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/5.Junction_permutation/2.logOR.ind.patient/100k_perm
echo "Log OR done" 


echo "Pseudobulking Data" 
#arg1 is filtered counts matrix 
#arg2 is filtered metadata
#arg3 is genotyping info
#arg4 is output from permutated observed difference 
#arg5 is output directory - combined matrix including pvalues, observed differences, annotations, and pseudobulked counts 
	#------------Pseudobulk counts and merge all data into one output ------# 
#Rscript ./scripts/pseudobulk_celltypes_counts.R /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/3.Counts_Filtered_Matrix/counts.filtered.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/4.Leafcutter_Annotation/P1_metadata_annotations.csv /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/5.Junction_permutation/100k_perm_v3/mut_wt_junction_permutation_table_R_output.Rdata /gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p1/6.Pseudobulk_Permutation_Output/observed_difference_100k_perm/
