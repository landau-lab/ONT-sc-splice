
##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##------------------------------- Generate a list of cryptic junctions to look into ------------------------------- ##

##Do this outside the script so that you have a fine with gene.name and junction
#junctions <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/4.Comb_patients_merge_counts/CH_combined/logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info.txt")
#junctions.sub <- junctions[grepl("GATA1|SLC25A39|ERGIC3|DPH5|TMEM141|DUSP12", junctions$gene),]
#rownames(junctions.sub) <- junctions.sub$intron_junction
#junctions.of.interest <- paste(junctions.sub$intron_junction)
#junctions.of.interest.final <- paste0(junctions.sub$gene,":",junctions.sub$intron_junction)

#load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/junctions.to.use.Rdata")
#load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/junctions.to.use.CH.full.Rdata")
load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/junctions.to.use.CH.5readThresh.Rdata")
print("junctions")
junctions.of.interest
print("junctions with gene names")
junctions.of.interest.final
print("gene list")
gene.list

##------------------------------- File paths and arguments ------------------------------- ##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ch.integrated.Rdata"
intron.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/CH259.CH305.merged.intron.metadata.txt"
junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/CH259.CH305.merged.count.matrix.txt"
genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/CH259.CH305.merged.genotype.info.with.cell.type.txt"

pseudotime.column <- "pseudotime_monocle3"
genotype.column <- "Genotype_1UMI"
step  <- 200
cluster.column = "CellType" 
cluster.types = c("HSPC","IMP","MEP","EP")

##------------------------------- load in seurat object, intron metadata junction count matrix and genotyping.info files ------------------------------- ##

load(seurat.object.path) #samples.integrated@meta.data

intron.meta <- read.table(intron.meta.file.path, sep="\t", header = T)
intron.meta$intron.junction <- paste0(intron.meta$chr,":",intron.meta$start,":",intron.meta$end,":",intron.meta$strand)
clusterIDs <- paste(intron.meta[which(intron.meta$intron.junction %in% junctions.of.interest),]$clusterID_5p)
intron.meta.sub <- intron.meta[which(intron.meta$clusterID_5p %in% clusterIDs),]
rownames(intron.meta.sub) <- intron.meta.sub$intron.junction

print("Loading count matrix")
#m <- read.table(junction.count.file.path, header = T)
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[-1]

print("Done loading count matrix")

##Remove this line when using the merged count matrix
#colnames(m) <- paste0(colnames(m),"_1")

genotyping.table <- read.table(genotyping.table.file.path)


##------------------------------- Load the metadata (subsetting for only genotyped cells) & remove NA's from pseudotime column ------------------------------- ##

meta = as.data.frame(samples.integrated@meta.data[which(samples.integrated@meta.data$CellType %in% cluster.types & samples.integrated@meta.data$Genotype_1UMI_10x != "NA"),])
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
m.sub = m.sub[,ord]

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
  intron.meta.sub.temp$dPSI <- (intron.meta.sub.temp$mut.psi - intron.meta.sub.temp$wt.psi)

  out = as.data.frame(intron.meta.sub.temp[,c("mut.bulk.counts","wt.bulk.counts","cluster.cov.mut","cluster.cov.wt","mut.psi","wt.psi","dPSI")])
  rownames(out) <- paste0(intron.meta.sub.temp$intron.junction)
  out = out[junctions.of.interest,]
  
  return(out)

})


#junctions.of.interest.final <- paste0(intron.meta.sub[junctions.of.interest,]$gene,":",intron.meta.sub[junctions.of.interest,]$intron.junction)
#dPSI.per.window = (do.call(cbind, win.dPSI))
#colnames(dPSI.per.window) = as.character(1:ncol(dPSI.per.window))
#out.zscores <- as.data.frame(t(apply(dPSI.per.window, 1, function(x) (x-mean(x))/sd(x))))


##------------------------------- Get ONT.expression per window for each gene involved ------------------------------- ##

load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ONT.seruat.objects/ch.ont.samples.integrated.rna.normscale.Rdata")

gene.list <- intersect(rownames(ch.ont.samples.integrated@assays$RNA@data), gene.list)

if (length(setdiff(gene.list, rownames(ch.ont.samples.integrated@assays$RNA@data))) >0 ){
  print("gene's not in Gene expression matrix - probably due to differences in annotation")
  setdiff(gene.list, rownames(ch.ont.samples.integrated@assays$RNA@data))
}

read.count.threshold <- 0
#read.count.threshold <- 1000
#read.count.threshold <- 1500

ont.expression.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  length(shared.cells)
  nrow(meta.temp)
 
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  wt <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="WT"),])
  mut <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="MUT"),])
  message(paste(length(wt),"WT cells"))
  message(paste(length(mut),"MUT cells"))
  message(x[1],":", x[2])

  ch.ont.samples.integrated@meta.data$window <- 0
  ch.ont.samples.integrated@meta.data[shared.cells,]$window <- 1
  ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data$window ==1 & ch.ont.samples.integrated@meta.data$nCount_RNA < read.count.threshold),]$window <- 0

  Idents(ch.ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ch.ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
  expression.values$RNA
  
  out = expression.values$RNA[,"1"]
  return(out)
  
  })
  
ont.expression.per.window <- as.data.frame(ont.expression.per.window)
rownames(ont.expression.per.window) <- gene.list
colnames(ont.expression.per.window) <- as.character(1:ncol(ont.expression.per.window))


ont.mut.expression.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  length(shared.cells)
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  mut <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="MUT"),])
  message(paste(length(mut),"MUT cells"))
  message(x[1],":", x[2])
  
  ch.ont.samples.integrated@meta.data$window <- 0
  ch.ont.samples.integrated@meta.data[mut,]$window <- 1
  ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data$window ==1 & ch.ont.samples.integrated@meta.data$nCount_RNA < read.count.threshold),]$window <- 0
  
  Idents(ch.ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ch.ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
  expression.values$RNA
  out = expression.values$RNA[,"1"]
  return(out)
  
})

ont.mut.expression.per.window <- as.data.frame(ont.mut.expression.per.window)
rownames(ont.mut.expression.per.window) <- gene.list
colnames(ont.mut.expression.per.window) <- as.character(1:ncol(ont.mut.expression.per.window))


ont.wt.expression.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  length(shared.cells)
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  wt <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="WT"),])
  message(paste(length(wt),"WT cells"))
  message(x[1],":", x[2])
  
  ch.ont.samples.integrated@meta.data$window <- 0
  ch.ont.samples.integrated@meta.data[wt,]$window <- 1
  ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data$window ==1 & ch.ont.samples.integrated@meta.data$nCount_RNA < read.count.threshold),]$window <- 0
  
  Idents(ch.ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ch.ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
  expression.values$RNA
  out = expression.values$RNA[,"1"]
  return(out)
  
})

ont.wt.expression.per.window <- as.data.frame(ont.wt.expression.per.window)
rownames(ont.wt.expression.per.window) <- gene.list
colnames(ont.wt.expression.per.window) <- as.character(1:ncol(ont.wt.expression.per.window))


##------------------------------- Get ONT raw read counts per window for each gene involved ------------------------------- ##

ont.mut.and.wt.raw.counts.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  wt <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="WT"),])
  mut <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="MUT"),])
  message(paste(length(wt),"WT cells"))
  message(paste(length(mut),"MUT cells"))
  message(paste(length(shared.cells),"total cells"))
  message(x[1],":", x[2])
  bulk.counts <- as.data.frame(rowSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[gene.list,shared.cells])))
  colnames(bulk.counts) <- c("counts")
  out = bulk.counts
  return(out)
  
})

ont.mut.and.wt.raw.counts.per.window <- as.data.frame(ont.mut.and.wt.raw.counts.per.window)
colnames(ont.mut.and.wt.raw.counts.per.window) <- as.character(1:ncol(ont.mut.and.wt.raw.counts.per.window))


ont.mut.raw.counts.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  mut <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="MUT"),])
  message(paste(length(mut),"MUT cells"))
  message(x[1],":", x[2])
  bulk.counts <- as.data.frame(rowSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[gene.list,mut])))
  colnames(bulk.counts) <- c("counts")
  out = bulk.counts
  return(out)
  
})

ont.mut.raw.counts.per.window <- as.data.frame(ont.mut.raw.counts.per.window)
colnames(ont.mut.raw.counts.per.window) <- as.character(1:ncol(ont.mut.raw.counts.per.window))

ont.wt.raw.counts.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  ch.ont.samples.integrated.meta.sub <- ch.ont.samples.integrated@meta.data[shared.cells,]
  wt <- rownames(ch.ont.samples.integrated.meta.sub[which(ch.ont.samples.integrated.meta.sub$Genotype_1WT =="WT"),])
  message(paste(length(wt),"WT cells"))
  message(x[1],":", x[2])
  bulk.counts <- as.data.frame(rowSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[gene.list,wt])))
  colnames(bulk.counts) <- c("counts")
  out = bulk.counts
  return(out)
  
})

ont.wt.raw.counts.per.window <- as.data.frame(ont.wt.raw.counts.per.window)
colnames(ont.wt.raw.counts.per.window) <- as.character(1:ncol(ont.wt.raw.counts.per.window))


ont.total.raw.counts.per.window = lapply(define.windows, function(x){

  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
  message(paste(length(shared.cells),"total cells"))
  message(x[1],":", x[2])
  bulk.counts <- as.data.frame(colSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[,shared.cells])))
  bulk.count <- sum(bulk.counts[,1])
  out = bulk.count
  return(out)
})

ont.total.raw.counts.per.window <- unlist(ont.total.raw.counts.per.window)


ont.raw.counts.per.cell.per.window.per.gene <- list()

for (gene in gene.list){
  
  ont.raw.counts.per.cell.per.window = lapply(define.windows, function(x){
    
    meta.temp = m.meta[c(x[1]:x[2]),]
    shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
    message(x[1],":", x[2])
    
    counts <- as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[gene,shared.cells])
    colnames(counts)[1] <- c("gene.counts")
    counts$total <- colSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[,shared.cells]))
    counts$nCount_RNA <- ch.ont.samples.integrated@meta.data[shared.cells,]$nCount_RNA
    counts$norm.counts <- (counts$gene.counts/counts$total) * 10000
    counts$ln.norm.counts <- log1p(counts$norm.counts)
    counts$seurat.norm.counts <- ch.ont.samples.integrated@assays$RNA@data[gene,shared.cells]
    counts$genotype <- ch.ont.samples.integrated@meta.data[shared.cells,]$Genotype_1WT
    
    counts$used <- 0
    counts[which(counts$nCount_RNA >= read.count.threshold),]$used <- 1
    
    out = counts
    return(out)
    
  })
  
  ont.raw.counts.per.cell.per.window.per.gene[[gene]] <- ont.raw.counts.per.cell.per.window
}



##------------------------------- Save all output ------------------------------- ##

#save(pseudotime.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.celltype.unique.bucket.100.Rdata")
save(pseudotime.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.mut.and.wt.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window, ont.total.raw.counts.per.window, ont.raw.counts.per.cell.per.window.per.gene, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.info.sig.in.any.5readThresh.bucket.200.nCountThresh0.Rdata")

#write.table(m[rownames(intron.meta.sub),],file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.count.matrix.sub.txt", quote = F, row.names = T, col.names = T)


