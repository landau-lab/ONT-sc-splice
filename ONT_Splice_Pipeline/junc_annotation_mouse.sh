#!/bin/bash

#SBATCH --job-name=sicelore
#SBATCH --mem=500G
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderror_%j.log

source /etc/profile.d/modules.sh
module load java/1.9
module load samtools

module load racon

PATH=$PATH:/gpfs/commons/home/ahawkins/spoa/build/bin/
PATH=$PATH:/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2
echo $PATH

##### input arguments/ requirements 
## working directory must contain a folder labeled input_files 
## within input_files folder there are 3 folders
## 1.ONT_fastq -starts out empty but unzipped fastq is copied here 
## 2.short_read_files - contains barcodes.tsv (must be unzipped)  and possorted_genome_bam.bam and possorted_genome_bam.bam.bai from cellranger outs folder

workdir=$1 #working directory where all output files will be 
fastq=$2 #full path to fastq
fastqdir=$3
sample_name=$4
run_files=$5 ##location where all scripts are located
#genotype_info=$6
#pattern=$7
#nperm=$8
##barcodes=$8 #full path to barcodes file 

echo $sample_name
cd $workdir

junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline_mouse.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin))
