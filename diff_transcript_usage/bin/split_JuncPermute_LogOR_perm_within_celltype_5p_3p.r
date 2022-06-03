library(tidyverse)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-s", "--split"),
    type = "character",
    help = "path to folder with split counts matrix"
  ),
  make_option(
    opt_str = c("-g", "--genotype_file"),
    type = "character",
    help = "path to genotype tsv"
  ),
  make_option(
    opt_str = c("-n", "--num_perm"),
    type = "integer",
    default - 100000,
    help = "Number of permutations to perform"
  ),
  make_option(
    opt_str = c("-p", "--pattern"),
    type = "character",
    help = "pattern used to identify cells in integrated data (ie. _1, _2)"
  ),
  make_option(
    opt_str = c("-f", "--output_file"),
    type = "character",
    help = "output file name"
  ),
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character",
    help = "path to output directory"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))


path.to.split = opt$split
path.to.genotype = opt$genotype_file
nperm = opt$num_perm
pattern=opt$pattern
output.dir = opt$output_file
output.file = opt$output_dir 

nperm = as.numeric(nperm)

## test data
# path.to.genotype = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/2.Genotype_info/mds_p1.genotype.info.with.cell.type.txt"
# path.to.split = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/diff_transcript_output/split_cluster_files/split_1"


path.to.three.matrix = paste(path.to.split, "/three_prime/counts_files", sep = "")
path.to.five.matrix = paste(path.to.split, "/five_prime/counts_files", sep= "")
path.to.three.data = paste(path.to.split, "/three_prime/data_tables", sep = "")
path.to.five.data = paste(path.to.split, "/five_prime/data_tables", sep = "")

#load counts data
setwd(path.to.three.matrix)
file = list.files(path.to.three.matrix)
three.mtx = read.table(file, stringsAsFactors = FALSE, header = TRUE)
three.mtx = as.matrix(three.mtx)

setwd(path.to.five.matrix)
file = list.files(path.to.five.matrix)
five.mtx = read.table(file, stringsAsFactors = FALSE, header = TRUE)
five.mtx = as.matrix(five.mtx)
message("Matrix loaded")

## load genotype table 
geno = as.data.frame(read.table(path.to.genotype))
genotype = geno[,c("Genotype_1UMI", "Cell.Assignment")]
rownames(genotype) = gsub(rownames(genotype),pattern = pattern,replacement = "")
sum(rownames(genotype) %in% colnames(three.mtx))
genotype = genotype[rownames(genotype) %in% colnames(three.mtx),]
genotype = genotype[which(genotype$Genotype_1UMI != "AMB"),]
genotype$Cell.Assignment <- lapply(genotype$Cell.Assignment, as.character)
genotype = genotype[!is.na(genotype$Cell.Assignment),]

message("Genotype loaded")

#get list of wt and mut cells 
wt = rownames(genotype)[genotype$Genotype_1UMI == "WT"]
mut = rownames(genotype)[genotype$Genotype_1UMI == "MUT"]

#filter data for reads >= 5, clusters with junctions > 1 based on total observations 
setwd(path.to.three.data)
file = list.files(path.to.three.data)
three.data = read.csv(file)

setwd(path.to.five.data)
file = list.files(path.to.five.data)
five.data = read.csv(file)

#get list of unique clusters with n > 1
three.prime.clusters = as.character(unique(three.data$five_prime_ID))
five.prime.clusters = as.character(unique(five.data$three_prime_ID))

message("Data loaded")

# calculate observed odds ratio 

three.obs.ratio = list()
five.obs.ratio = list()

for(clust in three.prime.clusters){
  subset = three.data[which(three.data$five_prime_ID == clust),]
  subset$mut.all = sum(subset$obs.mut)
  subset$wt.all = sum(subset$obs.wt)
  subset$mut.remain = subset$mut.all - subset$obs.mut
  subset$wt.remain = subset$wt.all - subset$obs.wt
  rownames(subset) = subset$intron_junction
  for (junc in rownames(subset)){
    three.obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
  }
}

for(clust in five.prime.clusters){
  subset = five.data[which(five.data$three_prime_ID == clust),]
  subset$mut.all = sum(subset$obs.mut)
  subset$wt.all = sum(subset$obs.wt)
  subset$mut.remain = subset$mut.all - subset$obs.mut
  subset$wt.remain = subset$wt.all - subset$obs.wt
  rownames(subset) = subset$intron_junction
  for (junc in rownames(subset)){
    five.obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
  }
}

three.obs.ratio.num = as.numeric(three.obs.ratio)
names(three.obs.ratio.num) = names(three.obs.ratio)

five.obs.ratio.num = as.numeric(five.obs.ratio)
names(five.obs.ratio.num) = names(five.obs.ratio)

message("Observed difference calculated")

three.temp = rep(0,length(three.obs.ratio))
five.temp = rep(0,length(five.obs.ratio))# Create an empty vector to update in each iteration

#initiate shuffled data frame, do this only once 
three.shf.data = data.frame(intron_junction = three.data$intron_junction, five_prime_ID = three.data$five_prime_ID)
five.shf.data = data.frame(intron_junction = five.data$intron_junction, three_prime_ID = five.data$three_prime_ID)

cell.types = as.character(unique(genotype$Cell.Assignment))


for(x in 0:nperm){
  
  set.seed(x)
  
  # permute within each cell type 
  shuffle.genotype = list()
  for(type in cell.types){
    geno_vec = as.character(genotype[which(genotype$Cell.Assignment == type),"Genotype_1UMI"])
    names(geno_vec) = rownames(genotype[which(genotype$Cell.Assignment == type),])
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
  
  three.shf.data$shf.wt = rowSums(three.mtx[,colnames(three.mtx) %in% wt])
  three.shf.data$shf.mut = rowSums(three.mtx[,colnames(three.mtx) %in% mut])
  three.shf.data$total.reads.per.junction = three.shf.data$shf.wt+three.shf.data$shf.mut
  
  five.shf.data$shf.wt = rowSums(five.mtx[,colnames(five.mtx) %in% wt])
  five.shf.data$shf.mut = rowSums(five.mtx[,colnames(five.mtx) %in% mut])
  five.shf.data$total.reads.per.junction = five.shf.data$shf.wt+five.shf.data$shf.mut
  
  three.shf.ratio = list()
  five.shf.ratio = list()
  
  #calculate shuffled odds ratio 
  for(clust in three.prime.clusters){
    subset = three.shf.data[which(three.shf.data$five_prime_ID == clust),]
    subset$mut.all = sum(subset$shf.mut)
    subset$wt.all = sum(subset$shf.wt)
    subset$mut.remain = subset$mut.all - subset$shf.mut
    subset$wt.remain = subset$wt.all - subset$shf.wt
    rownames(subset) = subset$intron_junction
    for (junc in subset$intron_junction){
      three.shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
    }
  }
  
  for(clust in five.prime.clusters){
    subset = five.shf.data[which(five.shf.data$three_prime_ID == clust),]
    subset$mut.all = sum(subset$shf.mut)
    subset$wt.all = sum(subset$shf.wt)
    subset$mut.remain = subset$mut.all - subset$shf.mut
    subset$wt.remain = subset$wt.all - subset$shf.wt
    rownames(subset) = subset$intron_junction
    for (junc in subset$intron_junction){
      five.shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
    }
  }
  
  three.shf.ratio.num = as.numeric(three.shf.ratio)
  names(three.shf.ratio.num) = names(three.shf.ratio)
  ## calculate difference between observed ratio and shuffled ratio (TRUE FALSE vector)
  three.shf.diff = abs(three.obs.ratio.num) > abs(three.shf.ratio.num)
  
  ## add shuffled difference to temp vector, if TRUE will add 1 
  three.temp = three.temp + three.shf.diff
  
  five.shf.ratio.num = as.numeric(five.shf.ratio)
  names(five.shf.ratio.num) = names(five.shf.ratio)
  five.shf.diff = abs(five.obs.ratio.num) > abs(five.shf.ratio.num)
  
  five.temp = five.temp + five.shf.diff
  
  #progress(value = x, max.value = nperm, progress.bar = T)
  Sys.sleep(0.01)
  if(x == nperm) cat(" Permuted differences calculated")  
}

# compute pvals based on rank of observed difference in shuffled values 
three.pvals = 1 - three.temp/(nperm + 1)
five.pvals = 1 - five.temp/(nperm + 1)

message("Creating final data frames")
three.final = data.frame(pvalue = three.pvals, three.log.odds.ratio = three.obs.ratio.num, intron_junction = names(three.pvals))
three.final = left_join(three.data, three.final, by = "intron_junction")

three.cluster.cov = three.final %>% group_by(five_prime_ID) %>% summarise(three.mut.cluster.cov = sum(obs.mut), three.wt.cluster.cov = sum(obs.wt))
three.final = left_join(three.final, three.cluster.cov, by= "five_prime_ID")

five.final = data.frame(pvalue = five.pvals, five.log.odds.ratio = five.obs.ratio.num, intron_junction = names(five.pvals))
five.final = left_join(five.data, five.final, by = "intron_junction")

five.cluster.cov = five.final %>% group_by(three_prime_ID) %>% summarise(five.mut.cluster.cov = sum(obs.mut), five.wt.cluster.cov = sum(obs.wt))
five.final = left_join(five.final, five.cluster.cov, by= "three_prime_ID")

setwd(output.dir)
three.filename = paste("./alt_three_prime/", output.file, ".csv", sep = "")
five.filename = paste("./alt_five_prime/", output.file, ".csv", sep = "")
write.csv(three.final, file = three.filename, quote = FALSE, row.names = FALSE)
write.csv(five.final, file = five.filename, quote = FALSE, row.names = FALSE)
message("DONE!")
