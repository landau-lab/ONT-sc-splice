### inputs counts, metadata, and genotype information and subsets each by cluster ID 
library(Matrix)
library(parallel)
library(tidyverse)
library(matrixStats)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-c", "--counts"),
    type = "character",
    help = "path to counts matrix"
  ),
  make_option(
    opt_str = c("-g", "--genotype_file"),
    type = "character",
    help = "path to genotype tsv"
  ),
  make_option(
    opt_str = c("-d", "--metadata"),
    type = "character",
    help = "path to metadata tsv"
  ),
  make_option(
    opt_str = c("-p", "--pattern"),
    type = "character",
    help = "pattern used to identify cells in integrated data (ie. _1, _2)"
  ),
  make_option(
    opt_str = c("-r", "--min_reads"),
    type = "double",
    default = 5,
    help = "number of minimum reads for splice junction to cover"
  ),
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character",
    help = "path to output directory"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

path.to.matrix = opt$counts
path.to.genotype = opt$genotype_file
path.to.metadata = opt$metadata
pattern = opt$pattern
min_reads = opt$min_reads
output.dir = opt$output_dir


# # Test data
# path.to.matrix = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/leafcutter_outputs/MDS_P1_output/MDS_P1_perind_numbers.counts.txt"
# path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p1.genotype.info.with.cell.type.txt"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/strand_adjusted_outputs/strand_adjusted_metadata.csv"
# pattern = "_1"

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#merge by junctions in both metadata and counts matrix 
counts = read.table(path.to.matrix, header = TRUE, row.names = 1)
metadata = read.csv(path.to.metadata)
shared = intersect(rownames(counts), metadata$intron_junction)
mtx = counts[shared,]
metadata.filt = metadata[which(metadata$intron_junction %in% shared),]
message("Matrix loaded")


#load genotype information
geno = as.data.frame(read.table(path.to.genotype))
genotype = geno$Genotype_1UMI
genotype = as.character(genotype)
names(genotype) = rownames(geno)
names(genotype) = gsub(names(genotype), pattern = pattern,replacement = "") ##change pattern based on patient being input
sum(names(genotype) %in% colnames(mtx))
genotype = genotype[names(genotype) %in% colnames(mtx)]
genotype = genotype[genotype %in% c("WT", "MUT")]
message("Genotype loaded")

wt = names(genotype)[genotype == "WT"]
mut = names(genotype)[genotype == "MUT"]

intron_junction = data.frame(intron_junction = rownames(mtx))
intron_junction = intron_junction %>% separate(intron_junction, into = c("chr", "start", "end", "strand"), sep = ":")
intron_junction$five_prime = 'NA'
intron_junction$three_prime = 'NA'
intron_junction = intron_junction %>% mutate(five_prime = ifelse(strand=="+", start, end),
                                               three_prime = ifelse(strand == "+",end, start))
intron_junction$five_prime_ID = intron_junction %>% group_by(chr, five_prime, strand) %>% group_indices()
intron_junction$three_prime_ID = intron_junction %>% group_by(chr, three_prime, strand) %>% group_indices()
intron_junction$five_prime_ID = gsub("^", "clu_", intron_junction$five_prime_ID)
intron_junction$three_prime_ID = gsub("^", "clu_", intron_junction$three_prime_ID)
intron_junction$intron_junction = paste(intron_junction$chr, intron_junction$start, intron_junction$end, intron_junction$strand, sep = ":")


#filter data for reads >= 5, clusters with junctions > 1 based on total observations 
data = intron_junction
data$obs.wt = rowSums(mtx[,colnames(mtx) %in% wt])
data$obs.mut = rowSums(mtx[,colnames(mtx) %in% mut])
data$total.reads.per.junction = data$obs.wt+data$obs.mut
data = data[which(data$total.reads.per.junction >= min_reads),]

#get list of unique clusters with n > 1
alt.three.clusters = data %>% group_by(five_prime_ID) %>% tally()
alt.three.clust.list = alt.three.clusters %>% filter(n != 1) %>% select(five_prime_ID)
final.alt.three.clust = as.character(alt.three.clust.list$five_prime_ID)

alt.five.clusters = data %>% group_by(three_prime_ID) %>% tally()
alt.five.clust.list = alt.five.clusters %>% filter(n != 1) %>% select(three_prime_ID)
final.alt.five.clust = as.character(alt.five.clust.list$three_prime_ID)

## filter for three prime junctions 
data.filt.three = data[which(data$five_prime_ID %in% final.alt.three.clust),]
covered_junc = as.character(data.filt.three$intron_junction)
mtx.filt.three = mtx[covered_junc,]

## filter for five prime junctions 
data.filt.five = data[which(data$three_prime_ID %in% final.alt.five.clust),]
covered_junc = as.character(data.filt.five$intron_junction)
mtx.filt.five = mtx[covered_junc,]
message("Data filtered")

setwd(output.dir)

split_clusters_three = chunk(final.alt.three.clust, 10)
split_clusters_five = chunk(final.alt.five.clust, 10)

for (i in 1:length(split_clusters_three)) {
  data.split = data.filt.three[which(data.filt.three$five_prime_ID %in% split_clusters_three[[i]]),]
  split_junc = as.character(data.split$intron_junction)
  mtx.split = mtx.filt.three[split_junc,]
  #rownames(mtx.split) = split_junc
  workdir = paste("./split_", i, "/three_prime/data_tables/", sep = "")
  filename = paste("data.filt", i, "csv", sep = ".")
  filename=paste(workdir, filename, sep = "")
  write.csv(data.split, filename, quote = FALSE, row.names = FALSE)
  
  workdir = paste("./split_", i, "/three_prime/counts_files/", sep = "")
  filename = paste("mtx.filt", i, "txt", sep = ".")
  filename = paste(workdir, filename, sep="")
  write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
}

for (i in 1:length(split_clusters_five)) {
  data.split = data.filt.five[which(data.filt.five$three_prime_ID %in% split_clusters_five[[i]]),]
  split_junc = as.character(data.split$intron_junction)
  mtx.split = mtx.filt.five[split_junc,]
  #rownames(mtx.split) = split_junc
  workdir = paste("./split_", i, "/five_prime/data_tables/", sep = "")
  filename = paste("data.filt", i, "csv", sep = ".")
  filename=paste(workdir, filename, sep = "")
  write.csv(data.split, filename, quote = FALSE, row.names = FALSE)
  
  workdir = paste("./split_", i, "/five_prime/counts_files/", sep = "")
  filename = paste("mtx.filt", i, "txt", sep = ".")
  filename = paste(workdir, filename, sep="")
  write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
}

  
  