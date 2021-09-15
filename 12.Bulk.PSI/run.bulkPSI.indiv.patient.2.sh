#!/bin/bash
#SBATCH --job-name=bulk.PSI
#SBATCH --partition=bigmem,pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=200g
#SBATCH --output=err.and.out/MDS_P6_stdout_%j.log

module load R/3.5.1

patientID=MDS_P6

Rscript bulkPSI.indiv.patient.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/"$patientID"/output_files/leafcutter_outputs/"$patientID"_output/"$patientID"_perind_numbers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/"$patientID"/output_files/leafcutter_outputs/"$patientID"_output/"$patientID"_all.introns.info.w.primary.annotations.txt output/"$patientID".all.intron.metadata.with.bulkPSI.ONT.tsv
