###############################
#### this script generates the input matrices to do differential transcript usage ######
# input is full counts matrix and metadata
# output is strand adjusted counts matrix and metadata
# argument 1 is the full counts matrix from leafcutter
# argument 2 is the full annotated matrix from leafcutter (with additional annotations from Lloyd's annotation script)
# argument 3 is the directory to output the files 
###############################

library(tidyverse)
library(dplyr)

#load directories to input files 
args = commandArgs(TRUE)
count.file = args[1]
metadata.file = args[2]
output.dir = args[3]

#load in full counts matrix 
counts = read.table(count.file)
#counts = read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p2_align_smartseq_nofilter_perind_numers.counts.txt")
#counts = read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p1_align_smartseq_nofilter_perind_numers.counts.txt")
print("counts matrix loaded...")

#load in metadata 
all_introns_meta = read.csv(metadata.file)
print(dim(all_introns_meta))
#load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p1_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.Rdata")
#as.data.frame(load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/p2_v2_align_smartseq_nofilter_gtf_basic_all.introns.info.Rdata"))
print("Adjusting for strandedness...")
#adjust for strandedness
#all_introns_meta = all_introns_meta %>% dplyr::mutate(five_prime = ifelse(strand=="+", start, end),
                                               #three_prime = ifelse(strand=="+", end, start))


all_introns_meta$five_prime = "NA"
all_introns_meta$three_prime = "NA"
all_introns_meta$five_prime[which(all_introns_meta$strand == "+")] = "start"
all_introns_meta$five_prime[which(all_introns_meta$strand != "+")] = "end"
all_introns_meta$three_prime[which(all_introns_meta$strand =="+")] = "end"
all_introns_meta$three_prime[which(all_introns_meta$strand !="+")] = "start"
print("Taking strandedness into account")

#Taking strandedness into account 
all_introns_meta$fivep_diff <- NA 
all_introns_meta$threep_diff <- NA
all_introns_meta[which(all_introns_meta$strand == "+"), ]$fivep_diff <- all_introns_meta[which(all_introns_meta$strand == "+"), ]$start_diff
all_introns_meta[which(all_introns_meta$strand == "+"), ]$threep_diff <- all_introns_meta[which(all_introns_meta$strand == "+"), ]$end_diff
all_introns_meta[which(all_introns_meta$strand == "-"),]$fivep_diff <- (all_introns_meta[which(all_introns_meta$strand == "-"), ]$end_diff)*-1
all_introns_meta[which(all_introns_meta$strand == "-"),]$threep_diff <- (all_introns_meta[which(all_introns_meta$strand == "-"), ]$start_diff)*-1
print("Strand adjustment done")
print(dim(all_introns_meta))
#subset metadata by those that have strand information and filter counts information for only those that have strand information 

#all_introns_meta <- all_introns_meta %>% dplyr::mutate(intron_junction = paste(chr, start, end, sep = ":"))
all_introns_meta$intron_junction = 'NA'
all_introns_meta$intron_junction = paste(all_introns_meta$chr, all_introns_meta$start, all_introns_meta$end, sep = ":")
all_introns_meta = all_introns_meta %>% filter(strand != 'NA')
#counts = counts[which(unlist(lapply(strsplit(rownames(counts),split= ":clu"), "[",1)) %in% all_introns_meta$intron_junction),]


print("Writing output files...")
setwd(output.dir)
#write output matrices 
#write.table(counts,"/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_analysis_files/strand_adjusted_matrix/p1/P1_strand_adjusted_counts.txt", quote = FALSE)
#write.table(counts, "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_analysis_files/strand_adjusted_matrix/p2/P2_strand_adjusted_counts.txt", quote = FALSE)
#write_tsv(all_introns_meta,"/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_analysis_files/strand_adjusted_matrix/p1/P1_strand_adjusted_metadata.txt")
#write_tsv(all_introns_meta, "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_analysis_files/strand_adjusted_matrix/p2/P2_strand_adjusted_metadata.txt")
write.table(counts, "strand_adjusted_counts.txt")
write.table(all_introns_meta, "strand_adjusted_metadata.txt")
print("Done")
