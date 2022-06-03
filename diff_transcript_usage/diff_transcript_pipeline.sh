#!/bin/bash

#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2
#SBATCH --cpus-per-task=10

module load R/3.6.0


### arguments input for this file

# output_dir
# scripts_dir
# metadata
# counts
# genotype_info
# pattern
# sample_name
# nperm
# min_reads

# set defaults 
scripts_dir="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/ONT_Splice_Pipeline/diff_transcript_usage"
nperm=100000
min_reads=5

# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

###############################################################
######## Step 1: Run Strand Adjustment of Metadata ############
### Input: metadata output from leafcutter annotation 
### Input: path to output files/ working directory 
### Output: Srand adjusted metadata 

## format: Rscript strand_adjustment.R <metadata> <output directory> 

###############################################################

cd $output_dir 
mkdir strand_adjusted_metadata

Rscript "$scripts_dir"/bin/strand_adjustment.R \
  $metadata \
  "$output_dir"/strand_adjusted_metadata

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

Rscript "$scripts_dir"/bin/split_clusters_v2.R \
  --counts $counts \
  --genotype_file $genotype \
  --metadata "$output_dir"/strand_adjusted_metadata/strand_adjusted_metadata.csv \
  --pattern $pattern \
  --min_reads $min_reads \
  --output_dir "$output_dir"/diff_transcript_output/split_cluster_files

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
permute_jobids+=($(sbatch --job-name="$sample_name" "$scripts_dir"/bin/run_split_perm_within_celltype_5p_3p.sh \
  "$output_dir"/diff_transcript_output/split_cluster_files/split_"$i" \
  $genotype \
  $nperm \
  $pattern \
  "$output_dir"/diff_transcript_output/split_cluster_output \
  output_"$i" \
  "$scripts_dir"/bin))
done 

###########################################################
####### Step 4: Merge final output into one file and merge with all annotation information 
### Input: scripts_dir path 
### Input: outputs directory where all files are stored
### Input: strand adjusted metadata 
### Input: Final outfile 

mkdir merge_final_output

#sbatch "$scripts_dir"/bin/run_merge_output.sh "$scripts_dir"/bin "$output_dir"/diff_transcript_output/split_cluster_output "$output_dir"/strand_adjusted_metadata/strand_adjusted_metadata.csv "$output_dir"/diff_transcript_output/merge_final_output

merge=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/bin/run_merge_output.sh \
  "$scripts_dir"/bin \
  "$output_dir"/diff_transcript_output/split_cluster_output \
  "$output_dir"/strand_adjusted_metadata/strand_adjusted_metadata.csv \
  "$output_dir"/diff_transcript_output/merge_final_output))

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
permute_celltypes_jobids+=($(sbatch --job-name="$sample_name" "$scripts_dir"/bin/submit_split_celltype.sh \
  "$output_dir"/diff_transcript_output/split_cluster_files/split_"$i" \
  $genotype \
  $nperm \
  $pattern \
  "$output_dir"/diff_transcript_output/split_cluster_celltype_output \
  output_"$i" \
  "$scripts_dir"/bin))
done

##########################################################
##### Step 6: Merge output from DTU from cell types 
### input: scripts_dir path 
### input: outputs directory where all files are stored
### input: strand adjusted metadata
### input: genotype 
### input: Final output 

mkdir merge_final_celltype_output

merge_celltype=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/bin/run_merge_celltype_output.sh \
  "$scripts_dir"/bin \
  "$output_dir"/diff_transcript_output/split_cluster_celltype_output \
  "$output_dir"/strand_adjusted_metadata/strand_adjusted_metadata.csv \
  $genotype \
  "$output_dir"/diff_transcript_output/merge_final_celltype_output))

echo "Done" 

