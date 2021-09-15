library(Matrix)
library(parallel)
#install.packages("svMisc")
#library(svMisc)


args = commandArgs(TRUE)
path.to.matrix = args[1]
path.to.genotype = args[2]
#nperm = args[3]
output.dir = args[4]

nperm = 100000

# Test data
#path.to.matrix = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p2/3.Counts_Filtered_Matrix/counts.filtered.txt"
#path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p2.genotype.info.with.cell.type.txt"

runJunctionPermutation = function(path.to.matrix, path.to.genotype, nperm){
  
  mtx = read.table(path.to.matrix)
  mtx = mtx[,-1] #this removes the total junction reads column that is present from counts filtering
  mtx = as.matrix(mtx)
  #mtx = as(mtx,"sparseMatrix")
  message("Matrix loaded")
  
  geno = as.data.frame(read.table(path.to.genotype))
  genotype = geno$Genotype_1UMI
  genotype = as.character(genotype)
  names(genotype) = rownames(geno)
  names(genotype) = gsub(names(genotype), pattern = "_1",replacement = "") ##change pattern based on patient being input
  sum(names(genotype) %in% colnames(mtx))
  genotype = genotype[names(genotype) %in% colnames(mtx)]
  genotype = genotype[genotype %in% c("WT", "MUT")]
  message("Genotype loaded")
  
  wt = names(genotype)[genotype == "WT"]
  mut = names(genotype)[genotype == "MUT"]
  
  obs.mut = rowSums(mtx[,colnames(mtx) %in% mut])
  obs.wt = rowSums(mtx[,colnames(mtx) %in% wt])
  
  total.reads.per.junction = obs.wt+obs.mut
  obs.mut = obs.mut/sum(obs.mut)*1e6 # Counts per million reads (CPM)
  obs.wt = obs.wt/sum(obs.wt)*1e6
  
  obs.difference = obs.mut-obs.wt 
  FC = obs.mut/obs.wt
  message("Observed difference calculated")
  
  temp = rep(0,length(obs.difference)) # Create an empty vector to update in each iteration
  
  for(x in 0:nperm){
    
    set.seed(x)
    
    orig.names = names(genotype)
    shuffle.genotype = sample(genotype, size = length(genotype), replace = F)
    names(shuffle.genotype) = orig.names
    
    wt = names(shuffle.genotype)[shuffle.genotype == "WT"]
    mut = names(shuffle.genotype)[shuffle.genotype == "MUT"]
    
    shf.mut = rowSums(mtx[,colnames(mtx) %in% mut])
    shf.wt = rowSums(mtx[,colnames(mtx) %in% wt])
    
    shf.mut = shf.mut/sum(shf.mut)*1e6
    shf.wt = shf.wt/sum(shf.wt)*1e6
    
    shf.difference = shf.mut-shf.wt
    
    shf.difference = abs(obs.difference) > abs(shf.difference)
    
    temp = temp + shf.difference
    
    #progress(value = x, max.value = nperm, progress.bar = T)
    Sys.sleep(0.01)
    if(x == nperm) cat(" Permuted differences calculated")  
  }
  
  pvals = 1 - temp/(nperm + 1)
  adjusted.pvals = p.adjust(pvals, method = "fdr")
  final = data.frame(pvalue = pvals, FC = FC, obs.diff.cpm = obs.difference, total.reads = total.reads.per.junction)
  return(final)
  message("Done!")
}

out = runJunctionPermutation(path.to.matrix = path.to.matrix, path.to.genotype = path.to.genotype, nperm)
message("Done!")

setwd(output.dir)
save(out, file = "mut_wt_junction_permutation_table_R_output.Rdata")