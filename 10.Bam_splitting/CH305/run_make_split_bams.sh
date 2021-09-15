#!/bin/bash

#SBATCH --job-name=submit_split_bams
#SBATCH --mem=5g

## arg 1 = workdir
## arg 2  = genotype 
## arg 3 = celltypes 
## arg 4 = bam 

sbatch /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/celltype_bams/scripts/make_celltype_bams_pattern2.sh  /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/celltype_bams/CH305 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_305.genotype.info.with.cell.type.txt ""HSPC" "EP" "IMP" "MEP" "MkP" "NP"" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH305/output_files/sicelore_outputs/CH305_consensus.sorted.tags.GE.bam
