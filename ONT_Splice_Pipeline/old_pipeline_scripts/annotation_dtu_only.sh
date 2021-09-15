#!/bin/bash

#SBATCH --job-name=submit_annotation
#SBATCH --mem=32g
#SBATCH --partition=pe2


workdir=$1 #working directory where all output files will be 
sample_name=$2
run_files=$3 ##location where all scripts are located
genotype_info=$4
pattern=$5
nperm=$6

### this script is for use if you have a failure in sicelore and need to start after sicelore and just do annotation and DTU

cd $workdir

## step 5: Junction calling and annotation 
#sbatch --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin

#junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin))

## step 6: Differential transcript usage (with permutations within cell types)
sbatch "$run_files"/diff_transcript_usage_scripts/diff_transcript_pipeline.sh "$workdir"/output_files $run_files "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_all.introns.info.w.primary.annotations.txt "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  $genotype_info $pattern $sample_name $nperm


#diff_transcript=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/diff_transcript_usage_scripts/diff_transcript_pipeline.sh "$workdir"/output_files $run_files "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_all.introns.info.w.primary.annotations.txt "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  $genotype_info $pattern $sample_name $nperm))


