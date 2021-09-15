## calculate log odds ratio for all clusters with >=2 junctions per cluster and total reads in genotyped cells >= 5
## permute by shuffling genotype within each patient 
## calculate one odds ratio by pseudobulking across all patients 

library(Matrix)
library(tidyverse)
library(matrixStats)

args = commandArgs(TRUE)
path.to.split = args[1]
path.to.genotype = args[2]
nperm = args[3]
patient.names = args[4]
output.dir = args[5]
output.file = args[6]

print(args)

nperm = as.numeric(nperm)

# Test data
# path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/2.Genotype_info"
# patient.names = c("CH259", "CH305")
# pattern = c("_1", "_2")
# path.to.split = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output/split_cluster_files/split_1"


path.to.three.matrix = paste(path.to.split, "/three_prime/counts_files", sep = "")
path.to.five.matrix = paste(path.to.split, "/five_prime/counts_files", sep= "")
path.to.three.data = paste(path.to.split, "/three_prime/data_tables", sep = "")
path.to.five.data = paste(path.to.split, "/five_prime/data_tables", sep = "")


patient.names = unlist(strsplit(patient.names, split = ","))
pattern = list()
for(i in 1:length(patient.names)){
  pattern[[i]]= paste("_", i, sep = "")
}
  
pattern = unlist(pattern)

#read in matrix as named list of matrix
setwd(path.to.three.matrix)
files = list.files(path.to.three.matrix)
three.mtx.list = lapply(files, function(x) read.table(file=x, header = TRUE))
names(three.mtx.list) = patient.names

setwd(path.to.five.matrix)
files = list.files(path.to.five.matrix)
five.mtx.list = lapply(files, function(x) read.table(file=x, header = TRUE))
names(five.mtx.list) = patient.names
message("Matrix loaded")

#load genotype information
genotype.list = list()

setwd(path.to.genotype)
i = 1
for (file in list.files(path.to.genotype)){
  geno = as.data.frame(read.table(file))
  geno = geno[!is.na(geno$Cell.Assignment_2),]
  genotype = geno[,c("Genotype_1UMI", "Cell.Assignment_2")] 
  rownames(genotype) = gsub(rownames(genotype),pattern = pattern[i],replacement = "")
  sum(rownames(genotype) %in% colnames(three.mtx.list[[i]]))
  genotype = genotype[rownames(genotype) %in% colnames(three.mtx.list[[i]]),]
  genotype = genotype[which(genotype$Genotype_1UMI != "AMB"),]
  genotype$Cell.Assignment_2 <- lapply(genotype$Cell.Assignment_2, as.character)
  genotype.list[[patient.names[i]]] = genotype
  i = i+1
}
message("Genotype loaded")

setwd(path.to.three.data)
files = list.files(path.to.three.data)
three.data.comb = read.csv(files[1])

setwd(path.to.five.data)
files = list.files(path.to.five.data)
five.data.comb = read.csv(files[1])

message("All Data Loaded")

for (patient in patient.names){
  three.mtx.list[[patient]]$patient = patient
  five.mtx.list[[patient]]$patient = patient
}

three.mtx = bind_rows(three.mtx.list)
five.mtx = bind_rows(five.mtx.list)

three.prime.clusters = as.character(unique(three.data.comb$five_prime_ID))
five.prime.clusters = as.character(unique(five.data.comb$three_prime_ID))


three.obs.ratio = list()

for(cluster in three.prime.clusters){
  subset = as.data.frame(three.data.comb[which(three.data.comb$five_prime_ID == cluster),])
  subset$mut.all = sum(subset$obs.mut)
  subset$wt.all = sum(subset$obs.wt)
  subset$mut.remain = subset$mut.all - subset$obs.mut
  subset$wt.remain = subset$wt.all - subset$obs.wt
  rownames(subset) = paste(subset$intron_junction, subset$five_prime_ID, sep = ":")
  for (junc in rownames(subset)){
    three.obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
  }
}

three.obs.ratio.num = as.numeric(three.obs.ratio)
names(three.obs.ratio.num) = names(three.obs.ratio)

five.obs.ratio = list()
for (cluster in five.prime.clusters){
  subset = as.data.frame(five.data.comb[which(five.data.comb$three_prime_ID == cluster),])
  subset$mut.all = sum(subset$obs.mut)
  subset$wt.all = sum(subset$obs.wt)
  subset$mut.remain = subset$mut.all - subset$obs.mut
  subset$wt.remain = subset$wt.all - subset$obs.wt
  rownames(subset) = paste(subset$intron_junction, subset$three_prime_ID, sep = ":")
  for (junc in rownames(subset)){
    five.obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
  }
}

five.obs.ratio.num = as.numeric(five.obs.ratio)
names(five.obs.ratio.num) = names(five.obs.ratio)

three.data.comb$alt_three_prime_intron_junction = paste(three.data.comb$intron_junction, three.data.comb$five_prime_ID, sep = ":")
five.data.comb$alt_five_prime_intron_junction = paste(five.data.comb$intron_junction, five.data.comb$three_prime_ID, sep = ":")

#create two data frames with final output
three.patient.output = data.frame(three.obs.logOR.ratio = three.obs.ratio.num, alt_three_prime_intron_junction = names(three.obs.ratio))
three.patient.output = left_join(three.patient.output, three.data.comb, by = "alt_three_prime_intron_junction")
rownames(three.patient.output) = three.patient.output$alt_three_prime_intron_junction

five.patient.output = data.frame(five.obs.logOR.ratio = five.obs.ratio.num, alt_five_prime_intron_junction = names(five.obs.ratio))
five.patient.output = left_join(five.patient.output, five.data.comb, by = "alt_five_prime_intron_junction")
rownames(five.patient.output) = five.patient.output$alt_five_prime_intron_junction

message("Observed difference calculated")

#initialize to calculate whether or not obs OR > shf OR
temp.three = rep(0,length(three.obs.ratio.num)) # Create an empty vector to update in each iteration
temp.five = rep(0,length(five.obs.ratio.num))

#initiate shuffled data frame, do this only once
three.shf.data.list = list()
five.shf.data.list = list()
for (patient in patient.names){
  three.shf.data.list[[patient]] = data.frame(intron_junction = three.data.comb$intron_junction, five_prime_ID=three.data.comb$five_prime_ID)
  five.shf.data.list[[patient]] = data.frame(intron_junction = five.data.comb$intron_junction, three_prime_ID=five.data.comb$three_prime_ID)
}

#run through for every permutation
cell.types.list = list()
for (patient in patient.names){
  cell.types.list[[patient]] = as.character(unique(genotype.list[[patient]]$Cell.Assignment_2))
}

for(x in 0:nperm){

  set.seed(x)

  #cycle through each patient and shuffle genotypes before calculating OR
  #use only filtered data frames based on previous filtering done
  
  for (patient in patient.names) {


    shuffle.genotype = list()
    for(type in cell.types.list[[patient]]){
      geno_vec = as.character(genotype.list[[patient]][which(genotype.list[[patient]]$Cell.Assignment_2 == type),"Genotype_1UMI"])
      names(geno_vec) = rownames(genotype.list[[patient]][which(genotype.list[[patient]]$Cell.Assignment_2 == type),])
      orig.names = names(geno_vec)
      shuffle.genotype[[type]] = sample(geno_vec, size = length(geno_vec), replace = F)
      names(shuffle.genotype[[type]]) = orig.names
    }
    geno = flatten(shuffle.genotype)
    orig.names = names(geno)
    shf.genotype = as.character(geno)
    names(shf.genotype) = orig.names

    wt = names(shf.genotype)[shf.genotype == "WT"]
    mut = names(shf.genotype)[shf.genotype == "MUT"]

    three.shf.data = three.shf.data.list[[patient]]

    three.shf.data$shf.wt = rowSums(three.mtx.list[[patient]][,colnames(three.mtx.list[[patient]]) %in% wt])
    three.shf.data$shf.mut = rowSums(three.mtx.list[[patient]][,colnames(three.mtx.list[[patient]]) %in% mut])
    three.shf.data$total.reads.per.junction = three.shf.data$shf.wt+three.shf.data$shf.mut
    three.shf.data.list[[patient]] = three.shf.data
    
    five.shf.data = five.shf.data.list[[patient]]
    
    five.shf.data$shf.wt = rowSums(five.mtx.list[[patient]][,colnames(five.mtx.list[[patient]]) %in% wt])
    five.shf.data$shf.mut = rowSums(five.mtx.list[[patient]][,colnames(five.mtx.list[[patient]]) %in% mut])
    five.shf.data$total.reads.per.junction = five.shf.data$shf.wt+five.shf.data$shf.mut
    five.shf.data.list[[patient]] = five.shf.data
    
  }
  
  for (patient in patient.names){
    three.shf.data.list[[patient]]$patient = patient
    five.shf.data.list[[patient]]$patient = patient
  }
  
  three.shf.data = bind_rows(three.shf.data.list)
  five.shf.data = bind_rows(five.shf.data.list)
  
  three.shf.data.comb = three.shf.data %>% group_by(intron_junction, five_prime_ID) %>% summarise(shf.wt = sum(shf.wt), shf.mut = sum(shf.mut), total.reads.per.junction = sum(total.reads.per.junction))
  five.shf.data.comb = five.shf.data %>% group_by(intron_junction, three_prime_ID) %>% summarise(shf.wt = sum(shf.wt), shf.mut = sum(shf.mut), total.reads.per.junction = sum(total.reads.per.junction))
  
  

  three.shf.ratio = list()
  
  #log OR for each junction after shuffling genotypes within each patient
  for(clust in three.prime.clusters){
    subset = as.data.frame(three.shf.data.comb[which(three.shf.data.comb$five_prime_ID == clust),])
    subset$mut.all = sum(subset$shf.mut)
    subset$wt.all = sum(subset$shf.wt)
    subset$mut.remain = subset$mut.all - subset$shf.mut
    subset$wt.remain = subset$wt.all - subset$shf.wt
    rownames(subset) = paste(subset$intron_junction, subset$five_prime_ID, sep = ":")
    for (junc in subset$intron_junction){
      three.shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
    }
  }

  three.shf.ratio.num = as.numeric(three.shf.ratio)
  names(three.shf.ratio.num) = names(three.shf.ratio)
  
  ## calculate the same for five prime
  
  five.shf.ratio = list()
  for(clust in five.prime.clusters){
    subset = as.data.frame(five.shf.data.comb[which(five.shf.data.comb$three_prime_ID == clust),])
    subset$mut.all = sum(subset$shf.mut)
    subset$wt.all = sum(subset$shf.wt)
    subset$mut.remain = subset$mut.all - subset$shf.mut
    subset$wt.remain = subset$wt.all - subset$shf.wt
    rownames(subset) = paste(subset$intron_junction, subset$three_prime_ID, sep = ":")
    for (junc in subset$intron_junction){
      five.shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
    }
  }
  
  five.shf.ratio.num = as.numeric(five.shf.ratio)
  names(five.shf.ratio.num) = names(five.shf.ratio)

  #create output data frame with results for each patient
  three.shf.patient.output = data.frame(three.shf.logOR.ratio = three.shf.ratio.num, alt_three_prime_intron_junction = names(three.shf.ratio))
  rownames(three.shf.patient.output) = three.shf.patient.output$alt_three_prime_intron_junction
  
  five.shf.patient.output = data.frame(five.shf.logOR.ratio = five.shf.ratio.num, alt_five_prime_intron_junction = names(five.shf.ratio))
  rownames(five.shf.patient.output) = five.shf.patient.output$alt_five_prime_intron_junction


  three.shf.diff = abs(three.obs.ratio.num) > abs(three.shf.ratio.num)
  five.shf.diff = abs(five.obs.ratio.num) > abs(five.shf.ratio.num)

  temp.three = temp.three + three.shf.diff
  temp.five = temp.five + five.shf.diff

  #progress(value = x, max.value = nperm, progress.bar = T)
  Sys.sleep(0.01)
  if(x == nperm) cat(" Permuted differences calculated")
}

pvals.three = 1 - temp.three/(nperm + 1)
pvals.five = 1 - temp.five/(nperm + 1)

message("Creating final data frames")
final.three = data.frame(pvalue = pvals.three, three.obs.logOR.ratio = three.obs.ratio.num, intron_junction = three.patient.output$intron_junction)
final.three = left_join(final.three,three.patient.output)

three.cluster.cov = final.three %>% group_by(five_prime_ID) %>% summarise(three.mut.cluster.cov = sum(obs.mut), three.wt.cluster.cov = sum(obs.wt))
final.three = left_join(final.three, three.cluster.cov, by= "five_prime_ID")

final.five = data.frame(pvalue = pvals.five, five.obs.logOR.ratio = five.obs.ratio.num, intron_junction = five.patient.output$intron_junction)
final.five = left_join(final.five,five.patient.output)

five.cluster.cov = final.five %>% group_by(three_prime_ID) %>% summarise(five.mut.cluster.cov = sum(obs.mut), five.wt.cluster.cov = sum(obs.wt))
final.five = left_join(final.five, five.cluster.cov, by= "three_prime_ID")

message("Done with permutations!")

message("writing output")
setwd(output.dir)
three.filename = paste("./alt_three_prime/", output.file, ".csv", sep = "")
five.filename = paste("./alt_five_prime/", output.file, ".csv", sep = "")
write.csv(final.three, file = three.filename, quote = FALSE, row.names = FALSE)
write.csv(final.five, file = five.filename, quote = FALSE, row.names = FALSE)

#save(out, file = "MDS_combined_mut_wt_junction_permutation_logOR_R_output.Rdata")
message("Done!!")

