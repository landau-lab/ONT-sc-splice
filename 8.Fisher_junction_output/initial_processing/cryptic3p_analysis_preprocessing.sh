#!/bin/bash
#SBATCH --job-name=cryptic3p_analysis
#SBATCH --partition=bigmem
#SBATCH --mail-type=NONE
#SBATCH --mem=400gb
#SBATCH --output=ch305_cryptic3p_analysis_stdout_%j.log
#SBATCH --ntasks=10

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RenvForKnowlesPCA

#Rscript cryptic3p_analysis_preprocessing.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p1_align_smartseq_nofilter_perind_numers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p1_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.Rdata /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt -o mds_p1

#Rscript cryptic3p_analysis_preprocessing.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p2_align_smartseq_nofilter_perind_numers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p2_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.Rdata /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p2.genotype.info.with.cell.type.txt -o mds_p2

#Rscript cryptic3p_analysis_preprocessing.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p3_align_smartseq_nofilter_perind_numers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p3_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.Rdata /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p3.genotype.info.with.cell.type.txt -o mds_p3

#Rscript cryptic3p_analysis_preprocessing.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/CH259_nofilter_perind_numers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/CH259_nofilter_gtf_basic_all.introns.info.Rdata /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt -o ch259

Rscript cryptic3p_analysis_preprocessing.R /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/CH305_nofilter_perind_numers.counts.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/CH305_nofilter_gtf_basic_all.introns.info.Rdata /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_305.genotype.info.with.cell.type.txt -o ch305
