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
celltype = args[6]


# Test data
# path.to.matrix = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/1.Counts_matrix"
# path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined_2WT/5_read_threshold/2.Genotype_info"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/3.Annotated_metadata"
# patient.names = c("CH259", "CH305")
# pattern = c("_1", "_2")

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

patient.names = unlist(strsplit(patient.names, split = ","))
pattern = list()
for(i in 1:length(patient.names)){
  pattern[[i]]= paste("_", i, sep = "")
}

pattern = unlist(pattern)


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

#load genotype information
genotype.list = list()
setwd(path.to.genotype)
i = 1
for (file in list.files(path.to.genotype)){
  geno = as.data.frame(read.table(file))
  geno = geno[!is.na(geno$Cell.Assignment_2),]
  geno$total = geno$num.wt.calls_10x + geno$num.mut.calls_10x
  geno$mut_frac = geno$num.mut.calls/geno$total
  geno$Final_Genotype = NA
  geno[which(geno$num.wt.calls_10x >= 2 & geno$num.mut.calls_10x ==0), "Final_Genotype"] = "WT"
  geno[which(geno$num.mut.calls > 0),"Final_Genotype"] = "MUT"
  
  genotype = geno[,c("Cell.Assignment_2", "Final_Genotype")] 
  rownames(genotype) = gsub(rownames(genotype),pattern = pattern[i],replacement = "")
  sum(rownames(genotype) %in% colnames(mtx.list[[i]]))
  genotype = genotype[rownames(genotype) %in% colnames(mtx.list[[i]]),]
  genotype$Cell.Assignment_2 <- lapply(genotype$Cell.Assignment_2, as.character)
  genotype = genotype[!is.na(genotype$Final_Genotype),]
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

#filename = paste(output.dir, "/combined.data.txt", sep = "")
#write.table(combined, filename)
#dim(combined)

message("patient cluster ID's combined")


data.comb.list = list()
for (patient in patient.names){
  wt = rownames(genotype.list[[patient]])[genotype.list[[patient]]$Final_Genotype == "WT" & genotype.list[[patient]]$Cell.Assignment_2 == celltype]
  mut = rownames(genotype.list[[patient]])[genotype.list[[patient]]$Final_Genotype == "MUT" & genotype.list[[patient]]$Cell.Assignment_2 == celltype]
  
  
  data = combined[which(combined$patient == patient),]
  data$obs.wt = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% wt])
  data$obs.mut = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% mut])
  data$total.reads.per.junction = data$obs.wt+data$obs.mut
  data.comb.list[[patient]] = data
  data.comb.list[[patient]]$patient = patient 
}

data.comb = bind_rows(data.comb.list)
#dim(data.comb)
#write.table(data.comb, filename)

data.comb.counts = data.comb %>% group_by(intron_junction) %>% summarise(obs.wt = sum(obs.wt), obs.mut = sum(obs.mut), total.reads.per.junction = sum(total.reads.per.junction))
data.comb = data.comb %>% select( -obs.wt, -obs.mut, -total.reads.per.junction, -patient)
data.comb = distinct(data.comb)
data.comb = left_join(data.comb, data.comb.counts)
data.comb = data.comb[which(data.comb$total.reads.per.junction >= 5),]
dim(data.comb)

comb.metadata = distinct(Reduce(full_join, metadata.list) %>% select(intron_junction, startClass, endClass))
#dim(comb.metadata)

data.comb = inner_join(data.comb, comb.metadata)
#idx = which(data.comb$startClass == "" & data.comb$endClass == "")
#data.comb[-idx,]
#dim(data.comb)

clusters = data.comb %>% group_by(five_prime_ID) %>% tally()
filt.clusters = clusters %>% filter(n != 1)
clust.three = filt.clusters$five_prime_ID
alt.clust.three = data.comb[which(data.comb$endClass == "not_main_3_prime" | data.comb$startClass == "not_main_3_prime"), "five_prime_ID"]
alt.clust.three = intersect(alt.clust.three, clust.three)
length(alt.clust.three)

clusters = data.comb %>% group_by(three_prime_ID) %>% tally()
filt.clusters = clusters %>% filter(n != 1)
clust.five = filt.clusters$three_prime_ID
alt.clust.five = data.comb[which(data.comb$endClass == "not_main_5_prime" | data.comb$startClass == "not_main_5_prime"), "three_prime_ID"]
alt.clust.five = intersect(alt.clust.five, clust.five)
length(alt.clust.five)

## filter for three prime junctions
data.filt.three = data.comb[which(data.comb$five_prime_ID %in% alt.clust.three),]
covered_junc = as.character(data.filt.three$intron_junction)

three_prime_list = list()
for (patient in patient.names){
  three_prime_list[[patient]] = mtx.list[[patient]][covered_junc,]
}


## filter for five prime junctions 
data.filt.five = data.comb[which(data.comb$three_prime_ID %in% alt.clust.five),]
covered_junc = as.character(data.filt.five$intron_junction)

five_prime_list = list()
for (patient in patient.names){
  five_prime_list[[patient]] = mtx.list[[patient]][covered_junc,]
}


message("data filtered")


setwd(output.dir)
split_clusters_three = chunk(alt.clust.three, 100)
split_clusters_five = chunk(alt.clust.five, 100)
message("clusters split")


for (i in 1:length(split_clusters_three)){
  for (patient in patient.names) {
    data.split = data.filt.three[which(data.filt.three$five_prime_ID %in% split_clusters_three[[i]]),]
    split_junc = as.character(data.split$intron_junction)
    mtx.split = mtx.list[[patient]][split_junc,]
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
    data.split = data.filt.five[which(data.filt.five$three_prime_ID %in% split_clusters_five[[i]]),]
    split_junc = as.character(data.split$intron_junction)
    mtx.split = mtx.list[[patient]][split_junc,]
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