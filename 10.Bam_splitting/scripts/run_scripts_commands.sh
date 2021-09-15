#workdir=$1 - path to folder for output
#metadata=$2 - path to genotype metadata table
#bam=$3 - path to bam file

#sbatch path/to/script /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/* /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/* /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/* 


#MDS_P1 - Ally Ran


#MDS_P2 - used Genotype_1UMI
#sbatch make_celltype_bams_MDS_P2.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P2 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p2.genotype.info.with.cell.type.ONT_10x_merged.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P2_allreads_consensus.sorted.tags.GE.bam

#MDS_P3 - used Genotype_1UMI
#sbatch make_celltype_bams_MDS_P3.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P3 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p3.genotype.info.with.cell.type.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P3_allreads_consensus.sorted.tags.GE.bam

#CH259 - used Genotype_1UMI
#sbatch make_celltype_bams_CH259.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH259 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/CH259_consensus.sorted.tags.GE.bam

#Merging Genotype_2UMI cells
#sbatch merge_bams_CH259_2UMI.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH259 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt

#CH305 - Ally Ran

#Merging Genotype_2UMI cells
#sbatch merge_bams_CH305_2UMI.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/CH305 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_305.genotype.info.with.cell.type.txt

#MDS_P4 - NOT USING THIS BAM FILE


#MDS_P5_1 - used Genotype_1UMI_mutUMIFrac0.2 
#sbatch make_celltype_bams_MDS_P5_1.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_1 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p5_1.genotype.info.with.cell.type.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P5_1_consensus.sorted.tags.GE.bam

#MDS_P5_2 - used Genotype_1UMI_mutUMIFrac0.2 
#sbatch make_celltype_bams_MDS_P5_2.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P5_2 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p5_2.genotype.info.with.cell.type.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P5_2_consensus.sorted.tags.GE.bam

#MDS_P6 - used Genotype_1UMI_mutUMIFrac0.2 
#sbatch make_celltype_bams_MDS_P6.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/MDS_P6 /gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p6.genotype.info.with.cell.type.txt /gpfs/commons/groups/landau_lab/SF3B1_splice_project/2.Bam_files/2.ONT_sicelore_output/MDS_P6_consensus.sorted.tags.GE.bam 
