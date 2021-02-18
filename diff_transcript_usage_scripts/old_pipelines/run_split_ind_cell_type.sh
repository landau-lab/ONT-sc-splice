#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --mem=100g

workdir=$1
run_files=$2
metadata=$3
counts=$4
genotype=$5
pattern=$6
sample_name=$7
nperm=$8

cd $workdir

mkdir split_cluster_celltype_output
mkdir split_cluster_celltype_output/alt_three_prime
mkdir split_cluster_celltype_output/alt_five_prime
mkdir logs

celltype=()
for i in {1..1000}; do
celltype+=($(sbatch --job-name="$sample_name" "$run_files"/bin/submit_split_celltype.sh "$workdir"/split_cluster_files/split_"$i" $genotype $nperm $pattern "$workdir"/split_cluster_celltype_output output_"$i" "$run_files"/bin))
done

mkdir merge_final_celltype_output

merge_celltype=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/bin/run_merge_celltype_output.sh "$run_files"/bin "$workdir"/split_cluster_celltype_output $metadata $genotype "$workdir"/merge_final_celltype_output))

echo "Done" 
