## calculate log odds ratio for all clusters with >=2 junctions per cluster and total reads in genotyped cells >= 5
## permute by shuffling genotype  

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

#path.to.split = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output/split_cluster_files/split_1"

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
  geno = geno[!is.na(geno$Cell.Assignment),]
  geno$Cell.Assignment= as.character(geno$Cell.Assignment)
  geno[which(geno$Cell.Assignment == "EP_1"), "Cell.Assignment"] = "EP"
  geno[which(geno$Cell.Assignment == "EP_2"), "Cell.Assignment"] = "EP"
  genotype = geno[,c("Genotype_1UMI", "Cell.Assignment")] ## change Cell.Assignment column for MDS patients
  rownames(genotype) = gsub(rownames(genotype),pattern = pattern[i],replacement = "")
  sum(rownames(genotype) %in% colnames(three.mtx.list[[i]]))
  genotype = genotype[rownames(genotype) %in% colnames(three.mtx.list[[i]]),]
  genotype = genotype[which(genotype$Genotype_1UMI != "AMB"),]
  genotype$Cell.Assignment <- lapply(genotype$Cell.Assignment, as.character)
  genotype.list[[patient.names[i]]] = genotype
  i = i+1
}
message("Genotype loaded")

setwd(path.to.three.data)
files = list.files(path.to.three.data)
three.data.list = lapply(files, function(x) read.csv(file = x))
names(three.data.list) = patient.names

setwd(path.to.five.data)
files = list.files(path.to.five.data)
five.data.list = lapply(files, function(x) read.csv(file = x))
names(five.data.list) = patient.names

message("All Data Loaded")

three.prime.clusters = as.character(unique(three.data.list[[1]]$five_prime_ID))
five.prime.clusters = as.character(unique(five.data.list[[1]]$three_prime_ID))

three.patient.OR = list()
five.patient.OR = list()

for (patient in patient.names){
  wt = rownames(genotype.list[[patient]])[genotype.list[[patient]]$Genotype_1UMI == "WT"]
  mut = rownames(genotype.list[[patient]])[genotype.list[[patient]]$Genotype_1UMI == "MUT"]

  #load filtered data
  three.data.filt = three.data.list[[patient]]
  five.data.filt = five.data.list[[patient]]
  message("Data processed")

  three.obs.ratio = list()

  #get log odds ratio for each cluster
  ## do three_prime_ID first
  for(cluster in three.prime.clusters){
    subset = three.data.filt[which(three.data.filt$five_prime_ID == cluster),]
    subset$mut.all = sum(subset$obs.mut)
    subset$wt.all = sum(subset$obs.wt)
    subset$mut.remain = subset$mut.all - subset$obs.mut
    subset$wt.remain = subset$wt.all - subset$obs.wt
    rownames(subset) = paste(subset$intron_junction, subset$five_prime_ID, sep = ":")
    for (junc in rownames(subset)){
      three.obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
    }
  }
  #get total cells and obs ratio
  three.obs.ratio.num = as.numeric(three.obs.ratio)
  names(three.obs.ratio.num) = names(three.obs.ratio)

  five.obs.ratio = list()
  for (cluster in five.prime.clusters){
    subset = five.data.filt[which(five.data.filt$three_prime_ID == cluster),]
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

  three.data.filt$alt_three_prime_intron_junction = paste(three.data.filt$intron_junction, three.data.filt$five_prime_ID, sep = ":")
  five.data.filt$alt_five_prime_intron_junction = paste(five.data.filt$intron_junction, five.data.filt$three_prime_ID, sep = ":")

  #create two data frames with final output
  three.patient.output = data.frame(three.obs.logOR.ratio = three.obs.ratio.num, alt_three_prime_intron_junction = names(three.obs.ratio))
  three.patient.output = left_join(three.patient.output, three.data.filt, by = "alt_three_prime_intron_junction")
  rownames(three.patient.output) = three.patient.output$alt_three_prime_intron_junction
  colnames(three.patient.output) = lapply(colnames(three.patient.output), function(x) paste(patient,x,sep="_"))
  three.patient.OR[[patient]] = three.patient.output

  five.patient.output = data.frame(five.obs.logOR.ratio = five.obs.ratio.num, alt_five_prime_intron_junction = names(five.obs.ratio))
  five.patient.output = left_join(five.patient.output, five.data.filt, by = "alt_five_prime_intron_junction")
  rownames(five.patient.output) = five.patient.output$alt_five_prime_intron_junction
  colnames(five.patient.output) = lapply(colnames(five.patient.output), function(x) paste(patient,x,sep="_"))
  five.patient.OR[[patient]] = five.patient.output
}

#reduce one output from all patients into one dataframe from three dataframes
three.all = three.patient.OR %>% map( ~ .x %>%
                            rownames_to_column('intron_junction'))
for (patient in patient.names){
  three.all[[patient]]$intron_junction = sub(":clu.*", "", three.all[[patient]]$intron_junction)
}
three.all = three.all %>% purrr::reduce(full_join, by = "intron_junction")

five.all = five.patient.OR %>% map( ~ .x %>%
                                        rownames_to_column('intron_junction'))
for (patient in patient.names){
  five.all[[patient]]$intron_junction = sub(":clu.*", "", five.all[[patient]]$intron_junction)
}
five.all = five.all %>% purrr::reduce(full_join, by = "intron_junction")

#calculate mean observed difference
three.mat.names = lapply(patient.names, function(x) paste(x,"three.obs.logOR.ratio", sep = "_"))
three.obs.mat = three.all[,unlist(three.mat.names)]

five.mat.names = lapply(patient.names, function(x) paste(x,"five.obs.logOR.ratio", sep = "_"))
five.obs.mat = five.all[,unlist(five.mat.names)]

## weight by total junction coverage
three.total.junction.coverage = three.all[,grep("total.reads.per.junction", colnames(three.all))]
five.total.junction.coverage = five.all[,grep("total.reads.per.junction", colnames(five.all))]
three.weights = three.total.junction.coverage/rowSums(three.total.junction.coverage)
five.weights = five.total.junction.coverage/rowSums(five.total.junction.coverage)

i = 1
three.weight.obs.OR = list()
while (i <= dim(three.obs.mat)[1]){
  three.weight.obs.OR[[i]] = rowWeightedMeans(as.matrix(three.obs.mat[i,]), as.numeric(three.weights[i,]))
  i = i +1
}
three.weight.obs.OR = as.numeric(three.weight.obs.OR)
names(three.weight.obs.OR) = three.all$intron_junction

i = 1
five.weight.obs.OR = list()
while (i <= dim(five.obs.mat)[1]){
  five.weight.obs.OR[[i]] = rowWeightedMeans(as.matrix(five.obs.mat[i,]), as.numeric(five.weights[i,]))
  i = i +1
}
five.weight.obs.OR = as.numeric(five.weight.obs.OR)
names(five.weight.obs.OR) = five.all$intron_junction

message("Observed difference calculated")

#initialize to calculate whether or not obs OR > shf OR
temp.three = rep(0,length(three.weight.obs.OR)) # Create an empty vector to update in each iteration
temp.five = rep(0,length(five.weight.obs.OR))

#initiate shuffled data frame, do this only once
three.shf.data.list = list()
five.shf.data.list = list()
for (patient in patient.names){
  three.shf.data.list[[patient]] = data.frame(intron_junction = three.data.list[[patient]]$intron_junction, five_prime_ID=three.data.list[[patient]]$five_prime_ID)
  five.shf.data.list[[patient]] = data.frame(intron_junction = five.data.list[[patient]]$intron_junction, three_prime_ID=five.data.list[[patient]]$three_prime_ID)
}

#run through for every permutation
cell.types.list = list()
for (patient in patient.names){
  cell.types.list[[patient]] = as.character(unique(genotype.list[[patient]]$Cell.Assignment))
}

for(x in 0:nperm){

  set.seed(x)

  #initialize list for each patient to create list of results
  three.shf.patient.OR= list()
  five.shf.patient.OR= list()

  #cycle through each patient and shuffle genotypes before calculating OR
  #use only filtered data frames based on previous filtering done
  for (patient in patient.names) {


    shuffle.genotype = list()
    for(type in cell.types.list[[patient]]){
      geno_vec = as.character(genotype.list[[patient]][which(genotype.list[[patient]]$Cell.Assignment == type),"Genotype_1UMI"])
      names(geno_vec) = rownames(genotype.list[[patient]][which(genotype.list[[patient]]$Cell.Assignment == type),])
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

    three.shf.ratio = list()

    #log OR for each junction after shuffling genotypes within each patient
    for(clust in three.prime.clusters){
      subset = three.shf.data[which(three.shf.data$five_prime_ID == clust),]
      subset$mut.all = sum(subset$shf.mut)
      subset$wt.all = sum(subset$shf.wt)
      subset$mut.remain = subset$mut.all - subset$shf.mut
      subset$wt.remain = subset$wt.all - subset$shf.wt
      rownames(subset) = paste(subset$intron_junction, subset$five_prime_ID, sep = ":")
      for (junc in subset$intron_junction){
        three.shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
      }
    }
    #message("odds ratio calculated")

    three.shf.ratio.num = as.numeric(three.shf.ratio)
    names(three.shf.ratio.num) = names(three.shf.ratio)

    ## calculate the same for five prime
    five.shf.data = five.shf.data.list[[patient]]

    five.shf.data$shf.wt = rowSums(five.mtx.list[[patient]][,colnames(five.mtx.list[[patient]]) %in% wt])
    five.shf.data$shf.mut = rowSums(five.mtx.list[[patient]][,colnames(five.mtx.list[[patient]]) %in% mut])
    five.shf.data$total.reads.per.junction = five.shf.data$shf.wt+five.shf.data$shf.mut

    five.shf.ratio = list()
    for(clust in five.prime.clusters){
      subset = five.shf.data[which(five.shf.data$three_prime_ID == clust),]
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
    colnames(three.shf.patient.output) = lapply(colnames(three.shf.patient.output), function(x) paste(patient,x,sep="_"))
    three.shf.patient.OR[[patient]] = three.shf.patient.output

    five.shf.patient.output = data.frame(five.shf.logOR.ratio = five.shf.ratio.num, alt_five_prime_intron_junction = names(five.shf.ratio))
    rownames(five.shf.patient.output) = five.shf.patient.output$alt_five_prime_intron_junction
    colnames(five.shf.patient.output) = lapply(colnames(five.shf.patient.output), function(x) paste(patient,x,sep="_"))
    five.shf.patient.OR[[patient]] = five.shf.patient.output
  }

  #merge output into one dataframe
  three.shf.all = three.shf.patient.OR %>% map( ~ .x %>%
                                      rownames_to_column('intron_junction'))
  for (patient in patient.names){
    three.shf.all[[patient]]$intron_junction = sub(":clu.*", "", three.shf.all[[patient]]$intron_junction)
  }
  three.shf.all = three.shf.all %>% purrr::reduce(full_join, by = "intron_junction")

  five.shf.all = five.shf.patient.OR %>% map( ~ .x %>%
                                                  rownames_to_column('intron_junction'))
  for (patient in patient.names){
    five.shf.all[[patient]]$intron_junction = sub(":clu.*", "", five.shf.all[[patient]]$intron_junction)
  }
  five.shf.all = five.shf.all %>% purrr::reduce(full_join, by = "intron_junction")

  #calculate mean observed difference
  three.mat.names = lapply(patient.names, function(x) paste(x,"three.shf.logOR.ratio", sep = "_"))
  three.shf.mat = three.shf.all[,unlist(three.mat.names)]

  five.mat.names = lapply(patient.names, function(x) paste(x,"five.shf.logOR.ratio", sep = "_"))
  five.shf.mat = five.shf.all[,unlist(five.mat.names)]
  
  ## weight by total junction coverage calculated from observed values 
  
  i = 1
  three.weight.shf.OR = list()
  while (i <= dim(three.shf.mat)[1]){
    three.weight.shf.OR[[i]] = rowWeightedMeans(as.matrix(three.shf.mat[i,]), as.numeric(three.weights[i,]))
    i = i +1
  }
  three.weight.shf.OR = as.numeric(three.weight.shf.OR)
  names(three.weight.shf.OR) = three.shf.all$intron_junction
  
  i = 1
  five.weight.shf.OR = list()
  while (i <= dim(five.shf.mat)[1]){
    five.weight.shf.OR[[i]] = rowWeightedMeans(as.matrix(five.shf.mat[i,]), as.numeric(five.weights[i,]))
    i = i +1
  }
  five.weight.shf.OR= as.numeric(five.weight.shf.OR)
  names(five.weight.shf.OR) = five.shf.all$intron_junction

  three.shf.diff = abs(three.weight.obs.OR) > abs(three.weight.shf.OR)
  five.shf.diff = abs(five.weight.obs.OR) > abs(five.weight.shf.OR)

  temp.three = temp.three + three.shf.diff
  temp.five = temp.five + five.shf.diff

  #progress(value = x, max.value = nperm, progress.bar = T)
  Sys.sleep(0.01)
  if(x == nperm) cat(" Permuted differences calculated")
}

pvals.three = 1 - temp.three/(nperm + 1)
pvals.five = 1 - temp.five/(nperm + 1)
final.three = data.frame(pvalue = pvals.three, three.weight.logOR.ratio = three.weight.obs.OR, intron_junction = three.all$intron_junction)
final.three = left_join(final.three,three.all, by = "intron_junction")
final.five = data.frame(pvalue = pvals.five, five.weight.logOR.ratio = five.weight.obs.OR, intron_junction = five.all$intron_junction)
final.five = left_join(final.five,five.all, by = "intron_junction")

message("Done with permutations!")

## add total cells with junction to output
# count.names = lapply(patient.names, function(x) paste(x,"total.cells", sep = "_"))
# cell.counts = out %>% select(unlist(count.names))
# for (i in 1:length(patient.names)) {
#   name = paste("MDS_P",i,sep = "")
#   cell.counts[,name] = ifelse(cell.counts[,i] == 0,0,1)
# }
# out$patients.with.junction = rowSums(cell.counts[, patient.names])

message("writing output")
setwd(output.dir)
three.filename = paste("./alt_three_prime/", output.file, ".csv", sep = "")
five.filename = paste("./alt_five_prime/", output.file, ".csv", sep = "")
write.csv(final.three, file = three.filename, quote = FALSE, row.names = FALSE)
write.csv(final.five, file = five.filename, quote = FALSE, row.names = FALSE)

#save(out, file = "MDS_combined_mut_wt_junction_permutation_logOR_R_output.Rdata")
message("Done!!")

