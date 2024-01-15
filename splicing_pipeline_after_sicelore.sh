#!/usr/bin/env bash
#SBATCH --job-name=splice_pipeline
#SBATCH --partition=pe2
#SBATCH --output=%x%j.log
#SBATCH --error=%x%j.log

######################################################################################
# This script is used to submit the GoT-ONT splice pipeline for an individual sample 

# The output will include differential junction usage between all MUT and WT cells 
# and between all MUT and WT cells across all specified cell types 

# Any of the following parameters can be invoked at the command line: 

# fastq : Full path to FASTQ file from ONT
# short_read_files : Full path to folder with short read files
# sample_name : Name of sample being processed 
# scripts_dir : Path to where the GoT ONT pipeline scripts are stored 
# genotype_info : Table with genotype calls for each cell. 
# pattern : final characters from a barcode to trim, example if the barcode is 
#    AGGGCTACAGTTG_2, the pattern is "_2"
# output_dir : Where to store all output files  
# nperm : Number of permutations to perform for permutation testing. Default is 100000. 
# min_reads : Minimum number of reads a junction must have across all cells to be 
#   included in differential transcript usage. Default is 5.
# sicelore_dir : Directory where sicelore is located. 
# minimap_dir : Directory where minimap2 is located. 

######################################################################################

############ Command line arguments and load modules #################################

# set defaults 
scripts_dir="/gpfs/commons/groups/landau_lab/mariela/tools/ONT-sc-splice/"
nperm=100000
min_reads=5

sample_name="MDS02A-chr21"
genotype_info="/gpfs/commons/groups/landau_lab/mariela/general_scripts/minimal_example/mds_p5_1.genotype.info.with.cell.type.txt" 
output_dir="/gpfs/commons/groups/landau_lab/mariela/general_scripts/minimal_example"


# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
        echo "${v}: $2"
    fi
    shift
done


############ Error Checks #########################################################
# build path for sicelore outputs 
sicelore_outputs="$output_dir/output_files/sicelore_outputs"

junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/junction_annotation/junc_calling_pipeline.w.exon.skip.sh \
  --output_location "${output_dir}/output_files" \
  --input_bam "${sicelore_outputs}/${sample_name}_consensus.sorted.tags.GE.bam" \
  --sample ${sample_name} \
  --scripts_dir "${scripts_dir}/junction_annotation/" ))

############## Diff Transcript Usage ####################################################

## step 6: Differential transcript usage (with permutations within cell types)

diff_transcript=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/diff_transcript_usage/diff_transcript_pipeline.sh \
  --output_dir "$output_dir"/output_files \
  --scripts_dir $scripts_dir/diff_transcript_usage \
  --metadata "$output_dir"/output_files/leafcutter_outputs/"$sample_name"_all.introns.info.w.primaryAnnotations.exonSkipping.txt \
  --counts "$output_dir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  \
  --genotype_info $genotype_info \
  --pattern "_2" \
  --sample_name $sample_name \
  --nperm $nperm \
  --min_reads $min_reads ))
