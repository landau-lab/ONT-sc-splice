#!/bin/bash

#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2
#SBATCH --cpus-per-task=10

module load R/3.6.0


### arguments input for this file

workdir=$1
run_files=$2
metadata=$3
counts=$4
genotype=$5
pattern=$6
sample_name=$7
nperm=$8


###############################################################
######## Step 1: Run Strand Adjustment of Metadata ############
### Input: metadata output from leafcutter annotation 
### Input: path to output files/ working directory 
### Output: Srand adjusted metadata 

## format: Rscript strand_adjustment.R <metadata> <output directory> 

###############################################################

cd $workdir 
mkdir strand_adjusted_metadata

Rscript "$run_files"/bin/strand_adjustment.R $metadata "$workdir"/strand_adjusted_metadata

##############################################################
######### Step 2: Split clusters for differential transcript usage
### Input: full counts matrix from leafcutter junction calling
### Input: Genotype matrix
### Input: Strand adjusted metadata
### Input: Pattern (i.e. barcode pattern, "_1", "_2", or "_3")
### Input: Output directory 

## format: Rscript split_clusters_v2.R <counts> <genotype> <metadata> <pattern> <output>

###############################################################

mkdir diff_transcript_output
cd ./diff_transcript_output

mkdir split_cluster_files
cd ./split_cluster_files

for i in {1..1000}
do 
mkdir split_"$i"
mkdir split_"$i"/three_prime
mkdir split_"$i"/five_prime
mkdir split_"$i"/three_prime/counts_files
mkdir split_"$i"/three_prime/data_tables
mkdir split_"$i"/five_prime/counts_files
mkdir split_"$i"/five_prime/data_tables
done

Rscript "$run_files"/bin/split_clusters_v2.R $counts $genotype "$workdir"/strand_adjusted_metadata/strand_adjusted_metadata.csv $pattern "$workdir"/diff_transcript_output/split_cluster_files

###########################################################
######## Step 3: Batch submit each split cluster for differential analysis
### Input: path to split files 
### Input: path to genotype matrix 
### Input: Number of permutations
### Input: output directory 
### Input: output file name 
### output: differential transcript table for 3p and 5p in two separate folders

## format: sbatch run_split_perm_within_celltype_5p_3p.sh <path to split files> <genotype> <nperm> <output.dir> <output.file> 

###########################################################

cd ..
mkdir split_cluster_output
mkdir split_cluster_output/alt_three_prime
mkdir split_cluster_output/alt_five_prime
mkdir logs

permute_jobids=()
for i in {1..1000}; do
permute_jobids+=($(sbatch --job-name="$sample_name" "$run_files"/bin/run_split_perm_within_celltype_5p_3p.sh "$workdir"/diff_transcript_output/split_cluster_files/split_"$i" $genotype $nperm $pattern "$workdir"/diff_transcript_output/split_cluster_output output_"$i" "$run_files"/bin))
done 

###########################################################
####### Step 4: Merge final output into one file and merge with all annotation information 
### Input: run_files path 
### Input: outputs directory where all files are stored
### Input: strand adjusted metadata 
### Input: Final outfile 

mkdir merge_final_output

#sbatch "$run_files"/bin/run_merge_output.sh "$run_files"/bin "$workdir"/diff_transcript_output/split_cluster_output "$workdir"/strand_adjusted_metadata/strand_adjusted_metadata.csv "$workdir"/diff_transcript_output/merge_final_output

merge=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/bin/run_merge_output.sh "$run_files"/bin "$workdir"/diff_transcript_output/split_cluster_output "$workdir"/strand_adjusted_metadata/strand_adjusted_metadata.csv "$workdir"/diff_transcript_output/merge_final_output))

##########################################################
####### Step 5: Run DTU for all cell types 
### input: path to split files (use the same split files from mut vs. wt) 
### input: path to genotype 
### input: number of perm 
### input: pattern 
### input: output directory for all files 
### input: run files 
### output: differential transcript table with pvalue and log(odds ratio) for each cell type in two separate folders

## format: sbatch submit_split_celltype.sh <path to split files> <genotype> <nperm> <pattern> <output.dir> <output.file> <run files> 

###########################################################

mkdir split_cluster_celltype_output
mkdir split_cluster_celltype_output/alt_three_prime
mkdir split_cluster_celltype_output/alt_five_prime

permute_celltypes_jobids=()
for i in {1..1000}; do
permute_celltypes_jobids+=($(sbatch --job-name="$sample_name" "$run_files"/bin/submit_split_celltype.sh "$workdir"/diff_transcript_output/split_cluster_files/split_"$i" $genotype $nperm $pattern "$workdir"/diff_transcript_output/split_cluster_celltype_output output_"$i" "$run_files"/bin))
done

##########################################################
##### Step 6: Merge output from DTU from cell types 
### input: run_files path 
### input: outputs directory where all files are stored
### input: strand adjusted metadata
### input: genotype 
### input: Final output 

mkdir merge_final_celltype_output

merge_celltype=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/bin/run_merge_celltype_output.sh "$run_files"/bin "$workdir"/diff_transcript_output/split_cluster_celltype_output "$workdir"/strand_adjusted_metadata/strand_adjusted_metadata.csv $genotype "$workdir"/diff_transcript_output/merge_final_celltype_output))

echo "Done" 

