##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##------------------------------- File paths and arguments ------------------------------- ##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ch.integrated.Rdata"

pseudotime.column <- "pseudotime_monocle3"
genotype.column <- "Genotype_1UMI"
step  <- 60
cluster.column = "CellType" 
cluster.types = c("HSPC","IMP","MEP","EP")


##------------------------------- load in objects ------------------------------- ##

print("Loading seurat object")
load(seurat.object.path) #samples.integrated@meta.data
print("Done loading seurat object")


##------------------------------- Grab metadata (subsetting for only genotyped cells) & remove NA's from pseudotime column ------------------------------- ##

meta = as.data.frame(samples.integrated@meta.data[which(samples.integrated@meta.data[,cluster.column] %in% cluster.types & samples.integrated@meta.data[,genotype.column] != "NA"),])
meta <- meta[!is.na(meta[,pseudotime.column]),]


##------------------------------- Order cells in pseudotime ------------------------------- ##

ord = rownames(meta)[order(meta[,pseudotime.column],decreasing = F)]
meta = meta[ord,]

##------------------------------- Get cells per sliding windows ------------------------------- ##
step.def = seq(from = 1, to = nrow(meta), by = step)
define.windows = list()

for(i in 1:(round(nrow(meta)/(step))-1)){
  define.windows[[i]] = c(step.def[i],step.def[i+10])
}

ind = lapply(define.windows, function(x) sum(is.na(x)) > 0)
define.windows = define.windows[!unlist(ind)]
#define.windows

##------------------------------- Get mean pseudotime value and cell type fractions per window ------------------------------- ##

pseudotime.mean = unlist(lapply(define.windows, function(x){
  meta.temp = m.meta[c(x[1]:x[2]),]
  mean(meta.temp[,pseudotime.column])
}))


ct.fractions <- list()
if(!is.null(cluster.column) & !is.null(cluster.types)){
  for (ct in cluster.types) {
    ct.fractions[[ct]] = unlist(lapply(define.windows, function(x){
      meta.temp = m.meta[c(x[1]:x[2]),]
      sum(meta.temp[,cluster.column] == ct) / nrow(meta.temp)
    }))
  }
}


##------------------------------- Get mut.fraction per window ------------------------------- ##

mut.fraction = unlist(lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(genotyping.table), rownames(meta.temp))
  genotyping.table.sub <- genotyping.table[shared.cells, ]
  wt <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="WT"),])
  mut <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT"),])
  out = length(mut)/(length(wt)+ length(mut))
  
}))


