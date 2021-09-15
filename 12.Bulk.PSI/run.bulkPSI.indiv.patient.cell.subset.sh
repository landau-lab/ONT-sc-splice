#!/bin/bash
#SBATCH --job-name=subset.bulk.PSI
#SBATCH --partition=bigmem,pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=200g
#SBATCH --output=err.and.out/CH502_DNMT3A.WT.only.stdout_%j.log

module load R/3.5.1

patientID=CH502_DNMT3A 

Rscript bulkPSI.indiv.patient.cell.subset.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/"$patientID"/output_files/leafcutter_outputs/"$patientID"_output/"$patientID"_perind_numbers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH502_DNMT3A/output_files/leafcutter_outputs/CH502.metadata.txt "Genotype_Approx_Match_No_Gene" "WT" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/"$patientID"/output_files/leafcutter_outputs/"$patientID"_output/"$patientID"_all.introns.info.w.primary.annotations.txt output/"$patientID".WT.only.all.intron.metadata.with.bulkPSI.ONT.tsv
