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
celltypes=$9 #format is ("HSPC" "MEP") 

echo $celltypes

for type in ${celltypes[*]};
do
echo $type

grep -v "$type" $genotype > "$workdir"/diff_transcript_output_remove_celltypes/"$type"/"$type"_genotype_table.txt
        celltype_genotype="$workdir"/diff_transcript_output_remove_celltypes/"$type"/"$type"_genotype_table.txt
done
