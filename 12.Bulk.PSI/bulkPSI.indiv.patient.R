#!/usr/bin/env Rscript
#options(echo=TRUE)
### ONT - Generating bulk PSI scores

##--------------------------- Load in necessary libraries ---------------------------##

library(tibble)
library(tidyverse)
library(optparse)

##--------------------------- Allow for the input files to be arguments and save them as variables ---------------------------##

option_parser=OptionParser(
  usage="%prog [options] counts.file.Rdata intron.metadata.file.txt output.file"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 3)

#counts.file <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_perind_numbers.counts.txt"
#intron.meta.data.file <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/strand_adjusted_metadata/strand_adjusted_metadata.csv"
#intron.meta.data.file <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_all.introns.info.w.primary.annotations.txt"

#counts.file <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/CH259.CH305.merged.count.matrix.txt"
#intron.meta.data.file <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/CH259.CH305.merged.intron.metadata.txt"


counts.file <- parsed_args$args[1]
intron.meta.data.file <- parsed_args$args[2]
output.file <- parsed_args$args[3]

##--------------------------- Merging the junction counts for patient across all cells  ---------------------------##


read.count.data <- read.table(counts.file, check.names=FALSE, header=T)
rownames(read.count.data) <- read.count.data[,1]
read.count.data <- read.count.data[,-1]

#intron.meta.data <- read.csv(intron.meta.data.file)

#For merged matricies
intron.meta.data <- read.csv(intron.meta.data.file, sep="\t")
intron.meta.data$X <- rownames(intron.meta.data)

intron.meta.data$intron.junction <- paste0(intron.meta.data$chr,":",intron.meta.data$start,":",intron.meta.data$end,":",intron.meta.data$strand)

##Subset to remove duplicated and only include shared junctions
list <- intron.meta.data[duplicated(intron.meta.data$intron.junction), ]$intron.junction
to_remove <- intron.meta.data[which(intron.meta.data$intron.junction %in% list & (intron.meta.data$strand != intron.meta.data$strand_1)),]$X

intron.meta.data <- intron.meta.data[which(!(intron.meta.data$X %in% to_remove)),]
rownames(intron.meta.data) <- intron.meta.data$intron.junction

read.count.data <- read.count.data[intersect(rownames(intron.meta.data), rownames(read.count.data)),]
intron.meta.data <- intron.meta.data[intersect(rownames(intron.meta.data), rownames(read.count.data)),]

##--------------------------- Calculating Bulk PSI values --------------------------- ##

intron.meta.data$total.junction.coverage <- rowSums(read.count.data)
intron.meta.data <- intron.meta.data %>% group_by(clusterID_5p) %>% mutate(cluster.coverage = sum(total.junction.coverage))
intron.meta.data <- intron.meta.data %>% group_by(clusterID_5p) %>% mutate(bulk.psi = (total.junction.coverage/cluster.coverage)*100)


##--------------------------- Save output table with bulk PSI scores ---------------------------##

write.table(intron.meta.data, file=output.file, quote = F, sep="\t")









