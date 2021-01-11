### inputs counts, metadata, and genotype information and subsets each by cluster ID 
library(Matrix)
library(parallel)
library(tidyverse)
library(matrixStats)
library(data.table)

args = commandArgs(TRUE)
path.to.matrix = args[1]
path.to.genotype = args[2]
path.to.metadata = args[3]
output.dir = args[4]
patient.names = args[5]
pattern = args[6]


# Test data
path.to.matrix = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/1.Counts_matrix"
path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/2.Genotype_info"
path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output/combined_metadata/combined_metadata.csv"
path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/3.Annotated_metadata"
patient.names = c("CH259", "CH305")
pattern = c("_1", "_2")

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

metadata = read.csv(path.to.metadata)

patient.names = unlist(strsplit(patient.names, split = ","))
pattern = unlist(strsplit(pattern, split = ","))


## load in matrix 
setwd(path.to.matrix)
files = list.files(path.to.matrix)
mtx.list = lapply(files, function(x) fread(file=x))
names(mtx.list) = patient.names
message("Matrix loaded")

setwd(path.to.metadata)
files = list.files(path.to.metadata)
metadata.list = lapply(files, function(x) read.csv(file=x))
names(metadata.list) = patient.names

for (patient in patient.names){
  mtx.list[[patient]] = setDF(mtx.list[[patient]])
  rownames(mtx.list[[patient]]) = mtx.list[[patient]]$intron_junction
  metadata.list[[patient]] = metadata.list[[patient]] %>% mutate(intron_junction = paste(chr, start, end, strand, sep = ":"))
  shared = intersect(rownames(mtx.list[[patient]]), metadata.list[[patient]]$intron_junction)
  mtx.list[[patient]] = mtx.list[[patient]][shared,]
  metadata.list[[patient]] = metadata.list[[patient]][which(metadata.list[[patient]]$intron_junction %in% shared),]
  print(dim(mtx.list[[patient]]))
  print(dim(metadata.list[[patient]]))
}

genotype.list = list()
#load genotype information
setwd(path.to.genotype)
i = 1
for (file in list.files(path.to.genotype)){
  geno = as.data.frame(read.table(file))
  geno = geno[!is.na(geno$Cell.Assignment_2),]
  genotype = geno$Genotype_1UMI
  genotype = as.character(genotype)
  names(genotype) = rownames(geno)
  names(genotype) = gsub(names(genotype), pattern = pattern[i],replacement = "") ##change pattern based on patient being input
  sum(names(genotype) %in% colnames(mtx.list[[i]]))
  genotype = genotype[names(genotype) %in% colnames(mtx.list[[i]])]
  genotype = genotype[genotype %in% c("WT", "MUT")]
  genotype.list[[patient.names[i]]] = genotype
  i = i+1
}
message("Genotype loaded")

##Create combined metadata with combined three prime ID and five prime ID







patient.counts.list = list()
for (patient in patient.names){
  ## for each patient add in three prime and five prime ID information to counts matrix and remove unstranded 
  
  intron_junctions = data.frame(intron_junction = rownames(mtx.list[[patient]]))
  intron_junctions = intron_junctions %>% separate(intron_junction, c("chr", "start", "end", "clusterID"), sep = ":") %>% select(chr, start, end)
  intron_junctions$start = as.numeric(intron_junctions$start)
  intron_junctions$end = as.numeric(intron_junctions$end)
  
  intron_junctions = left_join(intron_junctions, metadata[grep(patient, rownames(metadata)),], by = c("chr", "start", "end"))
  
  intron_junctions = intron_junctions %>% mutate(intron_junc = paste(chr, start, end, sep = ":")) %>%
                                                   mutate(three_prime_junctions = paste(chr,start, end, strand, three_prime_ID, sep = ":")) %>%
                                                   mutate(five_prime_junctions = paste(chr, start, end, strand, five_prime_ID, sep = ":")) %>%
                                                   mutate(unique_ID = paste(chr, start, end, strand, three_prime_ID, five_prime_ID, sep = ":"))
  intron_junctions = intron_junctions %>% filter(strand == "+" | strand == "-")
  
  patient.counts = mtx.list[[patient]]
  
  patient.counts$intron_junction = rownames(patient.counts)
  patient.counts = patient.counts %>% separate(intron_junction, into = c("chr", "start", "end", "clusterID"), sep = ":") %>%
    mutate(intron_junc = paste(chr, start, end, sep = ":"))
  shared = intersect(patient.counts$intron_junc, intron_junctions$intron_junc)
  patient.counts = patient.counts[which(patient.counts$intron_junc %in% shared),]
  patient.counts$start = as.numeric(patient.counts$start)
  patient.counts$end = as.numeric(patient.counts$end)
  patient.counts = left_join(patient.counts, intron_junctions, by = c("intron_junc", "chr", "start", "end"))
  patient.counts.list[[patient]] = patient.counts
}
message("patient cluster ID's combined")


## make one list of clusters to split for three prime IDs
clust.three = list()
for (patient in patient.names){
  clusters = metadata[which(metadata$patient == patient),]
  clusters = clusters %>% group_by(five_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.three[[patient]] = filt.clusters$five_prime_ID
}
clust.list.three = Reduce(intersect, clust.three)

clust.five = list()
for (patient in patient.names){
  clusters = metadata[which(metadata$patient == patient),]
  clusters = clusters %>% group_by(three_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.five[[patient]] = filt.clusters$three_prime_ID
}
clust.list.five = Reduce(intersect, clust.five)
message("Cluster lists aligned")

data.three.list = list()
data.five.list = list()
three_prime_list = list()
five_prime_list = list()

metadata = metadata %>% mutate(unique_ID = paste(chr, start, end, strand, three_prime_ID, five_prime_ID, sep = ":")) %>%
  mutate(intron_junction = paste(chr, start, end, strand, sep = ":"))
metadata$patient = as.character(metadata$patient)

for (patient in patient.names){
  wt = names(genotype.list[[patient]])[genotype.list[[patient]] == "WT"]
  mut = names(genotype.list[[patient]])[genotype.list[[patient]] == "MUT"]
  
  shared = intersect(rownames(mtx.list[[patient]]), metadata[which(metadata$patient==patient), "intron_junction"])
  
  #filter data for reads >= 5, clusters with junctions > 1 based on total observations 
  #data = data.frame(intron_junction = patient.counts.list[[patient]]$intron_junc, three_prime_ID = patient.counts.list[[patient]]$three_prime_ID, five_prime_ID = patient.counts.list[[patient]]$five_prime_ID,
                    #unique_ID = patient.counts.list[[patient]]$unique_ID)
  #data$obs.wt = rowSums(patient.counts.list[[patient]][,colnames(patient.counts.list[[patient]]) %in% wt])
  #data$obs.mut = rowSums(patient.counts.list[[patient]][,colnames(patient.counts.list[[patient]]) %in% mut])
  
  data = data.frame(intron_junction = shared, 
                    three_prime_ID = metadata[which(metadata$patient == patient & metadata$intron_junction %in% shared), "three_prime_ID"],
                    five_prime_ID = metadata[which(metadata$patient == patient & metadata$intron_junction %in% shared), "five_prime_ID"],
                    unique_ID = metadata[which(metadata$patient == patient & metadata$intron_junction %in% shared), "unique_ID"])
  shared = intersect(rownames(mtx.list[[patient]]), data$intron_junction)
  data[which(data$intron_junction %in% shared),]
  dim(mtx.list[[patient]][shared,])
  
  data$obs.wt = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% wt])
  data$obs.mut = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% mut])
  data$total.reads.per.junction = data$obs.wt+data$obs.mut
  data = data[which(data$total.reads.per.junction >= 5),]
  
  ## filter for three prime junctions
  data.filt.three = data[which(data$five_prime_ID %in% clust.list.three),]
  data.filt.three$patient = patient
  covered_junc = as.character(data.filt.three$unique_ID)
  three_prime_list[[patient]] = patient.counts.list[[patient]][which(patient.counts.list[[patient]]$unique_ID %in% covered_junc),]
  data.three.list[[patient]] = data.filt.three
  
  ## filter for five prime junctions 
  data.filt.five = data[which(data$three_prime_ID %in% clust.list.five),]
  data.filt.five$patient = patient
  covered_junc = as.character(data.filt.five$unique_ID)
  five_prime_list[[patient]] = patient.counts.list[[patient]][which(patient.counts.list[[patient]]$unique_ID %in% covered_junc),]
  data.five.list[[patient]] = data.filt.five
}
message("data filtered")

## create list of clusters found in at least 1 patient with total reads > 5 
all.data.three = do.call(rbind, data.three.list)
clust.three = list()
for (patient in patient.names) {
  clusters = all.data.three[grep(patient, rownames(all.data.three)),]
  clusters = clusters %>% group_by(five_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.three[[patient]] = filt.clusters$five_prime_ID
}

final.clust.list.three = Reduce(intersect, clust.three)

all.data.five = do.call(rbind, data.five.list)
clust.five = list()
for (patient in patient.names) {
  clusters = all.data.five[grep(patient, rownames(all.data.five)),]
  clusters = clusters %>% group_by(three_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.five[[patient]] = filt.clusters$three_prime_ID
}
final.clust.list.five = Reduce(intersect, clust.five)

setwd(output.dir)
split_clusters_three = chunk(final.clust.list.three, 250)
split_clusters_five = chunk(final.clust.list.five, 250)
messge("clusters split")

# test.mtx = list()
# test.data = list()
# for (patient in patient.names){
#   data.split = data.three.list[[patient]][which(data.three.list[[patient]]$three_prime_ID %in% split_clusters_three[[i]]),]
#   split_junc = as.character(data.split$intron_junction)
#   mtx.split = patient.counts.list[[patient]][which(patient.counts.list[[patient]]$intron_junc %in% split_junc),]
#   test.mtx[[patient]] = mtx.split
#   test.data[[patient]] = data.split
# }



for (i in 1:length(split_clusters_three)){
  for (patient in patient.names) {
    data.split = data.three.list[[patient]][which(data.three.list[[patient]]$five_prime_ID %in% split_clusters_three[[i]]),]
    split_junc = as.character(data.split$intron_junction)
    mtx.split = patient.counts.list[[patient]][which(patient.counts.list[[patient]]$intron_junc %in% split_junc),]
    #rownames(mtx.split) = split_junc
    workdir = paste("./split_", i, "/three_prime/data_tables/", sep = "")
    filename = paste(patient,"data.filt", i, "csv", sep = ".")
    filename=paste(workdir, filename, sep = "")
    write.csv(data.split, filename, quote = FALSE, row.names = FALSE)

    workdir = paste("./split_", i, "/three_prime/counts_files/", sep = "")
    filename = paste(patient,"mtx.filt", i, "txt", sep = ".")
    filename = paste(workdir, filename, sep="")
    write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
  } 
}


for (i in 1:length(split_clusters_five)){
  for (patient in patient.names) {
    data.split = data.five.list[[patient]][which(data.five.list[[patient]]$three_prime_ID %in% split_clusters_five[[i]]),]
    split_junc = as.character(data.split$intron_junction)
    mtx.split = patient.counts.list[[patient]][which(patient.counts.list[[patient]]$intron_junc %in% split_junc),]
    #rownames(mtx.split) = split_junc
    workdir = paste("./split_", i, "/five_prime/data_tables/", sep = "")
    filename = paste(patient,"data.filt", i, "csv", sep = ".")
    filename=paste(workdir, filename, sep = "")
    write.csv(data.split, filename, quote = FALSE, row.names = FALSE)

    workdir = paste("./split_", i, "/five_prime/counts_files/", sep = "")
    filename = paste(patient,"mtx.filt", i, "txt", sep = ".")
    filename = paste(workdir, filename, sep="")
    write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
  } 
}
message("Done splitting clusters!")


  
  