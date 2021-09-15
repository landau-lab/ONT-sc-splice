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
# path.to.matrix = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/1.Counts_matrix"
# path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/2.Genotype_info"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output/combined_metadata/combined_metadata.csv"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/3.Annotated_metadata"
# patient.names = c("CH259", "CH305")
# pattern = c("_1", "_2")

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

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
  geno = geno[!is.na(geno$Cell.Assignment),]
  genotype = geno$Genotype_1UMI_mutUMIFrac0.2
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
combined_data = list()
for (patient in patient.names){
  combined_data[[patient]] = data.frame(intron_junction = rownames(mtx.list[[patient]]))
  combined_data[[patient]] = combined_data[[patient]] %>% separate(intron_junction, into = c("chr", "start", "end", "strand_1"), sep = ":") %>%
    mutate(intron_junction = paste(chr, start, end, strand_1, sep = ":"))
}

combined = do.call(rbind, combined_data)
combined$patient = rownames(combined)
combined$patient = gsub("\\..*", "", combined$patient)

combined$five_prime = 'NA'
combined$three_prime = 'NA'
combined = combined %>% mutate(five_prime = ifelse(strand_1=="+", start, end),
                                             three_prime = ifelse(strand_1 == "+",end, start))
combined$five_prime_ID = combined %>% group_by(chr, five_prime, strand_1) %>% group_indices()
combined$three_prime_ID = combined %>% group_by(chr, three_prime, strand_1) %>% group_indices()
combined$five_prime_ID = gsub("^", "clu_", combined$five_prime_ID)
combined$three_prime_ID = gsub("^", "clu_", combined$three_prime_ID)

combined = combined %>% mutate(unique_ID = paste(chr, start, end, strand_1, three_prime_ID, five_prime_ID, sep = ":"))

message("patient cluster ID's combined")


data.three.list = list()
data.five.list = list()
three_prime_list = list()
five_prime_list = list()


for (patient in patient.names){
  wt = names(genotype.list[[patient]])[genotype.list[[patient]] == "WT"]
  mut = names(genotype.list[[patient]])[genotype.list[[patient]] == "MUT"]
  

  data = combined[which(combined$patient == patient),]
  data$obs.wt = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% wt])
  data$obs.mut = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% mut])
  data$total.reads.per.junction = data$obs.wt+data$obs.mut
  data = data[which(data$total.reads.per.junction >= 5),]
  
  clusters = data %>% group_by(five_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.three = filt.clusters$five_prime_ID
  
  clusters = data %>% group_by(three_prime_ID) %>% tally()
  filt.clusters = clusters %>% filter(n != 1)
  clust.five = filt.clusters$three_prime_ID
  
  ## filter for three prime junctions
  data.filt.three = data[which(data$five_prime_ID %in% clust.three),]
  covered_junc = as.character(data.filt.three$intron_junction)
  three_prime_list[[patient]] = mtx.list[[patient]][covered_junc,]
  data.three.list[[patient]] = data.filt.three
  
  ## filter for five prime junctions 
  data.filt.five = data[which(data$three_prime_ID %in% clust.five),]
  covered_junc = as.character(data.filt.five$intron_junction)
  five_prime_list[[patient]] = mtx.list[[patient]][covered_junc,]
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
message("clusters split")

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
    mtx.split = mtx.list[[patient]][split_junc,]
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
    mtx.split = mtx.list[[patient]][split_junc,]
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


  
  