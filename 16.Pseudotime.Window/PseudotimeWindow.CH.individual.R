
##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(Seurat)

##------------------------------- Generate a list of cryptic junctions to look into ------------------------------- ##

##Do this outside the script so that you have a fine with gene.name and junction
junctions <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/4.Comb_patients_merge_counts/CH_combined/logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info.txt")
junctions.sub <- junctions[grepl("GATA1|SLC25A39|ERGIC3|DPH5|TMEM141|DUSP12", junctions$gene),]
junctions.of.interest <- paste(junctions.sub$intron_junction)
junctions.of.interest.final <- paste0(junctions.sub$gene,":",junctions.sub$intron_junction)


##------------------------------- File paths and arguments ------------------------------- ##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ch.integrated.Rdata"
intron.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_all.introns.info.w.primary.annotations.txt"
junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_perind_numbers.counts.txt"
genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt"

pseudotime.column <- "pseudotime_monocle3"
genotype.column <- "Genotype_1UMI"
step  <- 50
cluster.column = "CellType" 
cluster.types = c("HSPC","IMP","MEP","EP")

##------------------------------- load in seurat object, intron metadata junction count matrix and genotyping.info files ------------------------------- ##

load(seurat.object.path) #samples.integrated@meta.data

intron.meta <- read.table(intron.meta.file.path, sep=",", header = T)
intron.meta$intron.junction <- paste0(intron.meta$chr,":",intron.meta$start,":",intron.meta$end,":",intron.meta$strand_1)
rownames(intron.meta) <- intron.meta$intron.junction
clusterIDs <- paste(intron.meta[which(intron.meta$intron.junction %in% junctions.of.interest),]$clusterID_5p)
intron.meta.sub <- intron.meta[which(intron.meta$clusterID_5p %in% clusterIDs),]

m <- read.table(junction.count.file.path, header = T)
rownames(m) <- m$intron_junction
m <- m[,-1]
##Remove this line when using the merged count matrix
colnames(m) <- paste0(colnames(m),"_1")

genotyping.table <- read.table(genotyping.table.file.path)


##------------------------------- Load the metadata (subsetting for only genotyped cells) & remove NA's from pseudotime column ------------------------------- ##

meta = as.data.frame(samples.integrated@meta.data[which(samples.integrated@meta.data$orig.ident == "CH_259" & samples.integrated@meta.data$CellType %in% cluster.types & samples.integrated@meta.data$Genotype_1UMI_10x != "NA"),])
meta <- meta[!is.na(meta[,pseudotime.column]),]

#Subset the metadata and junction count matrix to have the same cells (and remove any NAs)
shared.cells <- intersect(colnames(m), rownames(meta))
print("number of cell in long read data")
ncol(m)
print("number of cell in short read data")
nrow(meta)
print("number of shared cells")
length(shared.cells)

m.sub = m[,shared.cells]
m.meta = meta[shared.cells,]

ncol(m.sub)
nrow(m.meta)

##------------------------------- Order cells in pseudotime ------------------------------- ##

ord = rownames(m.meta)[order(m.meta[,pseudotime.column],decreasing = F)]
m.meta = m.meta[ord,]
m.sub = m[,ord]

# Get sliding windows
step.def = seq(from = 1, to = ncol(m.sub), by = step)
define.windows = list()

for(i in 1:(round(ncol(m.sub)/(step))-1)){
  define.windows[[i]] = c(step.def[i],step.def[i+10])
}

ind = lapply(define.windows, function(x) sum(is.na(x)) > 0)
define.windows = define.windows[!unlist(ind)]
define.windows


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


ct.fractions.2 <- as.data.frame(matrix(ncol = length(define.windows), nrow = length(cluster.types)))
colnames(ct.fractions.2) <- as.character(1:ncol(ct.fractions.2))
rownames(ct.fractions.2) <- cluster.types

if(!is.null(cluster.column) & !is.null(cluster.types)){
  for (ct in cluster.types) {
    ct.fractions.2[ct,] = unlist(lapply(define.windows, function(x){
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

##------------------------------- Get differential splicing score per window ------------------------------- ##

win.dPSI =  lapply(define.windows, function(x){
  
  temp = m.sub[,c(x[1]:x[2])]
  meta.temp = m.meta[c(x[1]:x[2]),]
  intron.meta.sub.temp = intron.meta.sub

  shared.cells <- intersect(rownames(genotyping.table), rownames(meta.temp))
  genotyping.table.sub <- genotyping.table[shared.cells, ]
  
  wt <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="WT"),])
  mut <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT"),])
  message(paste(length(wt),"WT cells"))
  message(paste(length(mut),"MUT cells"))
  message(x[1],":",x[2])
  
  intron.meta.sub.temp$mut.bulk.counts <- rowSums(temp[rownames(intron.meta.sub.temp),mut])
  intron.meta.sub.temp$wt.bulk.counts <- rowSums(temp[rownames(intron.meta.sub.temp),wt])
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(cluster.cov.mut=sum(mut.bulk.counts))
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(cluster.cov.wt=sum(wt.bulk.counts))
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(mut.psi = (mut.bulk.counts/cluster.cov.mut)*100)
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(wt.psi = (wt.bulk.counts/cluster.cov.wt)*100)
  
  intron.meta.sub.temp
  intron.meta.sub.temp$dPSI <- (intron.meta.sub.temp$mut.psi - intron.meta.sub.temp$wt.psi)
  out = as.data.frame(intron.meta.sub.temp[,"dPSI"])
  
  #out = as.data.frame(intron.meta.sub.temp[,c("mut.bulk.counts","wt.bulk.counts","cluster.cov.mut","cluster.cov.wt","mut.psi","wt.psi","dPSI")])
  rownames(out) <- paste0(intron.meta.sub.temp$gene,":",intron.meta.sub.temp$intron.junction)
  #out = out[junctions.of.interest.final,]
  return(out)

})


#junctions.of.interest.final <- paste0(intron.meta.sub[junctions.of.interest,]$gene,":",intron.meta.sub[junctions.of.interest,]$intron.junction)
dPSI.per.window = (do.call(cbind, win.dPSI))
colnames(dPSI.per.window) = as.character(1:ncol(dPSI.per.window))
#out.zscores <- as.data.frame(t(apply(dPSI.per.window, 1, function(x) (x-mean(x))/sd(x))))


##------------------------------- Get ONT.expression per window for each gene involved ------------------------------- ##

load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ONT.seruat.objects/ch.ont.samples.integrated.rna.normscale.Rdata")

gene.list <- paste(unique(intron.meta.sub$gene))

ont.expression.per.window = lapply(define.windows, function(x){
  
  #x <- unlist(define.windows[1])
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  length(shared.cells)
  nrow(meta.temp)
  wt <- rownames(ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data[shared.cells,"Genotype_1WT"] =="WT"),])
  mut <- rownames(ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data[shared.cells,"Genotype_1WT"] =="MUT"),])
  message(paste(length(wt),"WT cells"))
  message(paste(length(mut),"MUT cells"))
  
  message(x[1],":", x[2])

  ch.ont.samples.integrated@meta.data$window <- 0
  ch.ont.samples.integrated@meta.data[shared.cells,]$window <- 1

  Idents(ch.ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ch.ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
  expression.values$RNA
  
  out = expression.values$RNA[,"1"]
  return(out)
  
  })
  
ont.expression.per.window <- as.data.frame(ont.expression.per.window)
rownames(ont.expression.per.window) <- gene.list
colnames(ont.expression.per.window) <- as.character(1:ncol(ont.expression.per.window))


##------------------------------- Save all output ------------------------------- ##

save(pseudotime.mean, ct.fractions, ct.fractions.2, mut.fraction, dPSI.per.window, ont.expression.per.window, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/ch.259.info.Rdata")





