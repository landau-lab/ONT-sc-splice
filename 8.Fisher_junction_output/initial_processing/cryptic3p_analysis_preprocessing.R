require(tidyverse)
library(dplyr)
library(magrittr)
library(tibble)
library(optparse)
library(RANN)

option_parser=OptionParser(
  usage="%prog [options] <name>_perind_numers.counts.txt <name>.introns.info.Rdata <name>.genotype.info.with.cell.type.txt",
  option_list=list(
    make_option( c("-o","--output"), default="out", help="Path and prefix for output file"))
)

parsed_args <- parse_args(option_parser,  positional_arguments = 3)

counts_file <- parsed_args$args[1]
intron_meta_file <- parsed_args$args[2]
genotyping_file <- parsed_args$args[3]

output = parsed_args$options$output
results_file <- paste0(output ,"_all_introns_meta_file.Rdata" )

#Load in the intron_meta and counts files
counts_reads <- read.table(counts_file, check.names=FALSE)
load(intron_meta_file)

counts_reads_file <- counts_reads
all_introns_meta_file <- all_introns_meta
  
cell_names <- colnames(counts_reads_file)
num_cells <- ncol(counts_reads_file)

all_introns_meta_file = all_introns_meta_file %>% mutate(five_prime = ifelse(strand=="+", start, end),
                                                          three_prime = ifelse(strand=="+", end, start))

#Taking strandedness into account 
all_introns_meta_file$fivep_diff <- NA 
all_introns_meta_file$threep_diff <- NA

all_introns_meta_file[which(all_introns_meta_file$strand == "+"), ]$fivep_diff <- all_introns_meta_file[which(all_introns_meta_file$strand == "+"), ]$start_diff
all_introns_meta_file[which(all_introns_meta_file$strand == "+"), ]$threep_diff <- all_introns_meta_file[which(all_introns_meta_file$strand == "+"), ]$end_diff
  
all_introns_meta_file[which(all_introns_meta_file$strand == "-"),]$fivep_diff <- (all_introns_meta_file[which(all_introns_meta_file$strand == "-"), ]$end_diff)*-1
all_introns_meta_file[which(all_introns_meta_file$strand == "-"),]$threep_diff <- (all_introns_meta_file[which(all_introns_meta_file$strand == "-"), ]$start_diff)*-1

#Calculating the total coverage for each junction
all_introns_meta_file$Total_Cov <- rowSums(counts_reads_file)

##Filtering step for junctions with no strand ID and low coverage:
all_introns_meta_file <- all_introns_meta_file[which(!is.na(all_introns_meta_file$strand)),]
all_introns_meta_file <- all_introns_meta_file[which(all_introns_meta_file$Total_Cov >= 10),]

##Filter out junctions labelled as novel annotated pairs - these are exon skipped events
all_introns_meta_file <- all_introns_meta_file[which(all_introns_meta_file$verdict != "novel annotated pair"),]

##Filter out the junctions that have only 1 acceptor - we don't really care about those  
num_junc = all_introns_meta_file %>% group_by(gene, clusterID) %>% summarise(n = n()) %>% ungroup()
#plot(table(num_junc$n), xlab="# alt acceptor", ylab="count")
n_per_five_prime = all_introns_meta_file %>% left_join(num_junc, by=c("gene", "clusterID")) %$% n
all_introns_meta_file$cluSize <- n_per_five_prime
to_keep = n_per_five_prime > 1
all_introns_meta_file = all_introns_meta_file[to_keep,]

#num_junc_2 = all_introns_meta_file %>% group_by(gene, clusterID) %>% summarise(n = n()) %>% ungroup()
#plot(table(num_junc_2$n), xlab="# alt acceptor post filtering", ylab="count")

##Adding in the distances of 3P junctions between eachother (only recording the shortest distances)
all_introns_meta_file = all_introns_meta_file %>% group_by(gene, clusterID) %>% mutate( d = nn2(three_prime, k = 2)$nn.dists[,2] )

##Add the same rownames as in the counts_reads file
all_introns_meta_file <- as.data.frame(all_introns_meta_file)
rownames(all_introns_meta_file) <- paste(all_introns_meta_file$chr,":",all_introns_meta_file$start, ":",all_introns_meta_file$end,":",all_introns_meta_file$clusterID, sep="")

##Final Counts file to use:
counts_reads_file <- counts_reads_file[rownames(all_introns_meta_file ),]

##Differential Junction Usage (with PSI values) 

#Load in the full genotype information
full_genotyping_info <- read.table(genotyping_file)

rownames(full_genotyping_info) <- unlist(lapply(strsplit(rownames(full_genotyping_info),split= "_"), "[",1))
shared_cells <- intersect(rownames(full_genotyping_info), colnames(counts_reads_file))
genotyping_info_shared <- full_genotyping_info[shared_cells,]

genotype <- "Genotype_1UMI"

##Pseudobulk all MUT and WT cells
mutcells <- rownames(genotyping_info_shared[which(genotyping_info_shared[,genotype] =="MUT"),])
wtcells <- rownames(genotyping_info_shared[which(genotyping_info_shared[,genotype] =="WT"),])

all_introns_meta_file$mut_bulk_counts <- rowSums(counts_reads_file[,mutcells])
all_introns_meta_file$wt_bulk_counts <- rowSums(counts_reads_file[,wtcells])

#all_introns_meta_file <- all_introns_meta_file %>% group_by(gene, clusterID) %>% mutate(cluster_cov_all = sum(Total_Cov))
all_introns_meta_file <- all_introns_meta_file %>% group_by(gene, clusterID) %>% mutate(cluster_cov_mut_bulk = sum(mut_bulk_counts))
all_introns_meta_file <- all_introns_meta_file %>% group_by(gene, clusterID) %>% mutate(cluster_cov_wt_bulk = sum(wt_bulk_counts))

#Save the final file that's been filtered
save(all_introns_meta_file,version=2, file=results_file)

