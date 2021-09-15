#!/bin/bash
#SBATCH --job-name=bulk.PSI.merged
#SBATCH --partition=bigmem,pe2
#SBATCH --mail-type=NONE
#SBATCH --mem=200g
#SBATCH --output=err.and.out/MDS_P5.MDS_P6_stdout_%j.log

module load R/3.5.1


#Rscript bulkPSI.merged.patients.R "CH259","CH305" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/CH259.all.intron.metadata.with.bulkPSI.ONT.tsv /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/CH305.all.intron.metadata.with.bulkPSI.ONT.tsv output/CH259.CH305.all.intron.metadata.with.bulkPSI.ONT.tsv

#Rscript bulkPSI.merged.patients.R "MDS_P5_1","MDS_P5_2" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/MDS_P5_1.all.intron.metadata.with.bulkPSI.ONT.tsv /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/MDS_P5_2.all.intron.metadata.with.bulkPSI.ONT.tsv output/MDS_P5_1.MDS_P5_2.all.intron.metadata.with.bulkPSI.ONT.tsv

Rscript bulkPSI.merged.patients.R "MDS_P5_1.MDS_P5_2","MDS_P6" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/MDS_P5_1.MDS_P5_2.all.intron.metadata.with.bulkPSI.ONT.tsv /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/MDS_P6.all.intron.metadata.with.bulkPSI.ONT.tsv output/MDSP5.MDSP6.all.intron.metadata.with.bulkPSI.ONT.tsv

