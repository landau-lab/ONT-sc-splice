library(tidyverse)

args = commandArgs(TRUE)
path.to.split = args[1]
path.to.genotype = args[2]
nperm = args[3]
pattern=args[4]
output.dir = args[5]
output.file = args[6]

nperm = as.numeric(nperm)

## test data
#path.to.genotype = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/2.Genotype_info/mds_p1.genotype.info.with.cell.type.txt"
#path.to.split = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/diff_transcript_output/split_cluster_files/split_1"
#nperm=20

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


geno = as.data.frame(read.table(path.to.genotype))
geno$Cell.Assignment = as.character(geno$Cell.Assignment)
geno[grep("EP_", geno$Cell.Assignment), "Cell.Assignment"] = "EP"
genotype = geno[,c("Genotype_1UMI", "Cell.Assignment")]
rownames(genotype) = gsub(rownames(genotype),pattern = pattern,replacement = "")
sum(rownames(genotype) %in% colnames(three.mtx))
genotype = genotype[rownames(genotype) %in% colnames(three.mtx),]
genotype = genotype[which(genotype$Genotype_1UMI != "AMB"),]
#genotype$Cell.Assignment <- lapply(genotype$Cell.Assignment, as.character)
genotype[is.na(genotype)] <- "None"
message("Genotype loaded")


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


cell.types = as.character(unique(genotype$Cell.Assignment))

three.celltype.OR = three.data
five.celltype.OR = five.data

for (type in cell.types){
  
  #get list of wt and mut cells 
  wt = rownames(genotype)[(genotype$Genotype_1UMI == "WT" & genotype$Cell.Assignment == type)]
  mut = rownames(genotype)[(genotype$Genotype_1UMI == "MUT" & genotype$Cell.Assignment == type)]
  
  if (length(wt) >= 2 & length(mut) >= 2){
    
    ### calculate new obs mut and obs wt for each cell type
    three.data$obs.wt = rowSums(three.mtx[,colnames(three.mtx) %in% wt])
    three.data$obs.mut = rowSums(three.mtx[,colnames(three.mtx) %in% mut])
    three.data$total.reads.per.junction = three.data$obs.wt+three.data$obs.mut
    
    five.data$obs.wt = rowSums(five.mtx[,colnames(five.mtx) %in% wt])
    five.data$obs.mut = rowSums(five.mtx[,colnames(five.mtx) %in% mut])
    five.data$total.reads.per.junction = five.data$obs.wt+five.data$obs.mut
    
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
    
    #calculate total cluster coverage 
    three.cluster.cov = three.data %>% dplyr::group_by(five_prime_ID) %>% dplyr::summarise(three.mut.cluster.cov = sum(obs.mut), three.wt.cluster.cov = sum(obs.wt))
    three.merge = left_join(three.data, three.cluster.cov, by = "five_prime_ID")
    
    five.cluster.cov = five.data %>% dplyr::group_by(three_prime_ID) %>% dplyr::summarise(five.mut.cluster.cov = sum(obs.mut), five.wt.cluster.cov = sum(obs.wt))
    five.merge = left_join(five.data, five.cluster.cov, by = "three_prime_ID")
    
    #create data frame of final output
    three.patient.output = data.frame(three.obs.logOR.ratio = three.obs.ratio.num, obs.wt = three.data$obs.wt, obs.mut = three.data$obs.mut, 
                                      total.reads.per.junction = three.data$total.reads.per.junction, three.wt.cluster.cov = three.merge$three.wt.cluster.cov, 
                                      three.mut.cluster.cov = three.merge$three.mut.cluster.cov)
    colnames(three.patient.output) = lapply(colnames(three.patient.output), function(x) paste(type,x,sep="_"))
    three.patient.output$intron_junction = names(three.obs.ratio)
    three.celltype.OR = full_join(three.patient.output, three.celltype.OR)
    
    five.patient.output = data.frame(five.obs.logOR.ratio = five.obs.ratio.num,
                                     obs.wt = five.data$obs.wt, obs.mut = five.data$obs.mut, total.reads.per.junction = five.data$total.reads.per.junction, five.wt.cluster.cov = five.merge$five.wt.cluster.cov, 
                                     five.mut.cluster.cov = five.merge$five.mut.cluster.cov)
    colnames(five.patient.output) = lapply(colnames(five.patient.output), function(x) paste(type,x,sep="_"))
    five.patient.output$intron_junction = names(five.obs.ratio)
    five.celltype.OR = full_join(five.patient.output, five.celltype.OR)
  }
}

message("Observed difference calculated")

three.temp = list()
for (type in cell.types){
  three.temp[[type]] = rep(0,dim(three.celltype.OR)[1])
} # Create an empty vector to update in each iteration

five.temp = list()
for (type in cell.types){
  five.temp[[type]] = rep(0,dim(five.celltype.OR)[1])
}


#initiate shuffled data frame, do this only once 
three.shf.data = data.frame(intron_junction = three.data$intron_junction, five_prime_ID = three.data$five_prime_ID)
five.shf.data = data.frame(intron_junction = five.data$intron_junction, three_prime_ID = five.data$three_prime_ID)

for(x in 0:nperm){
  
  set.seed(x)
  
  for (type in cell.types){
    wt = rownames(genotype)[(genotype$Genotype_1UMI == "WT" & genotype$Cell.Assignment == type)]
    mut = rownames(genotype)[(genotype$Genotype_1UMI == "MUT" & genotype$Cell.Assignment == type)]
    
    if (length(wt) >=2 & length(mut) >=2){
      ## shuffle genotypes
      geno_vec = paste(genotype$Genotype_1UMI, genotype$Cell.Assignment, sep = "_")
      names(geno_vec) = rownames(genotype)
      orig.names = names(geno_vec)
      shuffle.genotype = sample(geno_vec, size = length(geno_vec), replace = F)
      names(shuffle.genotype) = orig.names
      
      wt.geno = paste("WT", type, sep = "_")
      mut.geno = paste("MUT", type, sep = "_")
      wt = names(shuffle.genotype)[shuffle.genotype == wt.geno]
      mut = names(shuffle.genotype)[shuffle.genotype == mut.geno]
      
      ## calculate shuffled wt/mut for each permutation 
      three.shf.data$shf.wt = rowSums(three.mtx[,colnames(three.mtx) %in% wt])
      three.shf.data$shf.mut = rowSums(three.mtx[,colnames(three.mtx) %in% mut])
      three.shf.data$total.reads.per.junction = three.shf.data$shf.wt+three.shf.data$shf.mut
      
      five.shf.data$shf.wt = rowSums(five.mtx[,colnames(five.mtx) %in% wt])
      five.shf.data$shf.mut = rowSums(five.mtx[,colnames(five.mtx) %in% mut])
      five.shf.data$total.reads.per.junction = five.shf.data$shf.wt+five.shf.data$shf.mut
      
      ## calculate log OR for each permutation 
      three.shf.ratio = list()
      five.shf.ratio = list()
      
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
      select = as.character(paste(type, "three.obs.logOR.ratio", sep = "_"))
      three.obs.ratio.num = as.numeric(three.celltype.OR[,colnames(three.celltype.OR) == select])
      three.shf.diff = abs(three.obs.ratio.num) > abs(three.shf.ratio.num)
      
      three.temp[[type]] = three.temp[[type]] + three.shf.diff
      
      five.shf.ratio.num = as.numeric(five.shf.ratio)
      names(five.shf.ratio.num) = names(five.shf.ratio)
      select = as.character(paste(type, "five.obs.logOR.ratio", sep = "_"))
      five.obs.ratio.num = as.numeric(five.celltype.OR[,colnames(five.celltype.OR) == select])
      five.shf.diff = abs(five.obs.ratio.num) > abs(five.shf.ratio.num)
      
      five.temp[[type]] = five.temp[[type]] + five.shf.diff
    }
  }
  Sys.sleep(0.01)
  if(x == nperm) cat(" Permuted differences calculated")
}

three.pvals = data.frame(intron_junction = three.celltype.OR$intron_junction)
five.pvals = data.frame(intron_junction = five.celltype.OR$intron_junction)
for (type in cell.types){
  name = paste(type, "pvalue", sep = "_")
  three.pvals[,name] = 1 - three.temp[[type]]/(nperm + 1)
  
  name = paste(type, "pvalue", sep = "_")
  five.pvals[,name] = 1 - five.temp[[type]]/(nperm + 1)
}

three.final = inner_join(three.celltype.OR, three.pvals)
five.final = inner_join(five.celltype.OR, five.pvals)

setwd(output.dir)
three.filename = paste("./alt_three_prime/", output.file, ".csv", sep = "")
five.filename = paste("./alt_five_prime/", output.file, ".csv", sep = "")
write.csv(three.final, file = three.filename, quote = FALSE, row.names = FALSE)
write.csv(five.final, file = five.filename, quote = FALSE, row.names = FALSE)
message("DONE!")
