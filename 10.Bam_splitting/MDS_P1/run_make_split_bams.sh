#!/bin/bash

#SBATCH --job-name=submit_split_bams
#SBATCH --mem=5g

## arg 1 = workdir
## arg 2  = genotype 
## arg 3 = celltypes 
## arg 4 = bam 

sbatch /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/celltype_bams/scripts/make_celltype_bams.sh /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/celltype_bams/MDS_P1 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt ""CC" "HSPC" "EP_1" "EP_2" "HSPC_EP" "IMP" "MEP" "MkP" "NP"" /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_combined/output_files/sicelore_outputs/MDS_P1_merged_consensus.sorted.tags.GE.bam
