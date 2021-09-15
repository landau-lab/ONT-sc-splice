
##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##------------------------------- Generate a list of cryptic junctions to look into ------------------------------- ##


load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/junctions.to.use.CH.dPSI2.Rdata")
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
step  <- 60
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

meta = as.data.frame(samples.integrated@meta.data[which(samples.integrated@meta.data$CellType %in% cluster.types & samples.integrated@meta.data$Genotype_1UMI_10x == "MUT"),])
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

#mut.fraction = unlist(lapply(define.windows, function(x){
  
#  meta.temp = m.meta[c(x[1]:x[2]),]
#  shared.cells <- intersect(rownames(genotyping.table), rownames(meta.temp))
#  genotyping.table.sub <- genotyping.table[shared.cells, ]
#  wt <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="WT"),])
#  mut <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT"),])
#  out = length(mut)/(length(wt)+ length(mut))
#  
#}))

##------------------------------- Get differential splicing score per window ------------------------------- ##

win.mutpsi =  lapply(define.windows, function(x){
  
  temp = m.sub[,c(x[1]:x[2])]
  meta.temp = m.meta[c(x[1]:x[2]),]
  intron.meta.sub.temp = intron.meta.sub

  shared.cells <- intersect(rownames(genotyping.table), rownames(meta.temp))
  message(paste(length(shared.cells),"All cells"))
  genotyping.table.sub <- genotyping.table[shared.cells, ]
  
  #wt <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="WT"),])
  mut <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT"),])
  #message(paste(length(wt),"WT cells"))
  message(paste(length(mut),"MUT cells"))
  message(x[1],":",x[2])
  
  intron.meta.sub.temp$mut.bulk.counts <- rowSums(temp[rownames(intron.meta.sub.temp),mut])
  #intron.meta.sub.temp$wt.bulk.counts <- rowSums(temp[rownames(intron.meta.sub.temp),wt])
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(cluster.cov.mut=sum(mut.bulk.counts))
  #intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(cluster.cov.wt=sum(wt.bulk.counts))
  intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(mut.psi = (mut.bulk.counts/cluster.cov.mut)*100)
  #intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(wt.psi = (wt.bulk.counts/cluster.cov.wt)*100)
  #intron.meta.sub.temp$dPSI <- (intron.meta.sub.temp$mut.psi - intron.meta.sub.temp$wt.psi)

  out = as.data.frame(intron.meta.sub.temp[,c("mut.bulk.counts","cluster.cov.mut","mut.psi")])
  rownames(out) <- paste0(intron.meta.sub.temp$intron.junction)
  out = out[junctions.of.interest,]
  
  return(out)

})



##------------------------------- Get per ct counts per window ------------------------------- ##


win.ct.counts =  lapply(define.windows, function(x){
  
  temp = m.sub[,c(x[1]:x[2])]
  meta.temp = m.meta[c(x[1]:x[2]),]
  intron.meta.sub.temp = intron.meta.sub
  
  shared.cells <- intersect(rownames(genotyping.table), rownames(meta.temp))
  message(paste(length(shared.cells),"All cells"))
  genotyping.table.sub <- genotyping.table[shared.cells, ]
  
  mut <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT"),])
  message(paste(length(mut),"MUT cells"))
  message(x[1],":",x[2])
  
  #intron.meta.sub.temp$mut.bulk.counts <- rowSums(temp[rownames(intron.meta.sub.temp),mut])
  #intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(cluster.cov.mut=sum(mut.bulk.counts))
  #intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(mut.psi = (mut.bulk.counts/cluster.cov.mut)*100)
  
  for (ct in cluster.types) {
    
    ct.mut.cells <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT" & genotyping.table.sub[,"Cell.Assignment_2"] == ct),])
    message(paste(length(ct.mut.cells), ct ,"MUT cells"))
    message(ct.mut.cells[1:5])
    
    if (length(ct.mut.cells) <= 1){
      message(paste("too few", ct ,"MUT cells"))
      intron.meta.sub.temp[,paste0(ct,".mut.bulk.counts")] <- NA
      intron.meta.sub.temp[,paste0(ct,".cluster.cov.mut")] <- NA
      intron.meta.sub.temp[,paste0(ct,".mut.psi")] <- NA
      
      
    } else {
      
      intron.meta.sub.temp[,paste0(ct,".mut.bulk.counts")] <- rowSums(temp[rownames(intron.meta.sub.temp),ct.mut.cells])
      
      if (ct == "HSPC"){
        message("Adding HSPC counts")
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(HSPC.cluster.cov.mut=sum(HSPC.mut.bulk.counts))
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(HSPC.mut.psi = (HSPC.mut.bulk.counts/HSPC.cluster.cov.mut)*100)
        rownames(intron.meta.sub.temp) <- rownames(intron.meta.sub)
        #print(intron.meta.sub.temp[130:150,c("HSPC.mut.bulk.counts","HSPC.cluster.cov.mut","HSPC.mut.psi")])
        
        next
      } 
      
      if (ct =="IMP"){
        message("Adding IMP counts")
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(IMP.cluster.cov.mut=sum(IMP.mut.bulk.counts))
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(IMP.mut.psi = (IMP.mut.bulk.counts/IMP.cluster.cov.mut)*100)
        rownames(intron.meta.sub.temp) <- rownames(intron.meta.sub)
        #print(intron.meta.sub.temp[130:150,c("IMP.mut.bulk.counts","IMP.cluster.cov.mut","IMP.mut.psi")])
        
        next
      }
      
      if (ct =="MEP"){
        message("Adding MEP counts")
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(MEP.cluster.cov.mut=sum(MEP.mut.bulk.counts))
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(MEP.mut.psi = (MEP.mut.bulk.counts/MEP.cluster.cov.mut)*100)
        rownames(intron.meta.sub.temp) <- rownames(intron.meta.sub)
        next
      }
      
      if (ct =="EP"){
        message("Adding EP counts")
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(EP.cluster.cov.mut=sum(EP.mut.bulk.counts))
        intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(EP.mut.psi = (EP.mut.bulk.counts/EP.cluster.cov.mut)*100)
        rownames(intron.meta.sub.temp) <- rownames(intron.meta.sub)
      }
      
    }
    
  }
  
  num.columns <- ncol(intron.meta.sub.temp)
  num.cols.to.keep <- length(cluster.types)*3
  out = as.data.frame(intron.meta.sub.temp[,(num.columns-num.cols.to.keep):num.columns])
  rownames(out) <- paste0(intron.meta.sub.temp$intron.junction)
  out = out[junctions.of.interest,]
  
  return(out)
  
})



##------------------------------- Get ONT.expression per window for each gene involved ------------------------------- ##

load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ONT.seruat.objects/ch.ont.samples.integrated.rna.normscale.Rdata")

gene.list <- intersect(rownames(ch.ont.samples.integrated@assays$RNA@data), gene.list)

if (length(setdiff(gene.list, rownames(ch.ont.samples.integrated@assays$RNA@data))) >0 ){
  print("gene's not in Gene expression matrix - probably due to differences in annotation")
  setdiff(gene.list, rownames(ch.ont.samples.integrated@assays$RNA@data))
}


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
  #ch.ont.samples.integrated@meta.data[which(ch.ont.samples.integrated@meta.data$window ==1),]$window <- 0
  
  Idents(ch.ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ch.ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
  expression.values$RNA
  out = expression.values$RNA[,"1"]
  return(out)
  
})

ont.mut.expression.per.window <- as.data.frame(ont.mut.expression.per.window)
rownames(ont.mut.expression.per.window) <- gene.list
colnames(ont.mut.expression.per.window) <- as.character(1:ncol(ont.mut.expression.per.window))



##------------------------------- Get ONT raw read counts per window for each gene involved ------------------------------- ##

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



#ont.raw.counts.per.cell.per.window.per.gene <- list()

#for (gene in gene.list){
#  
#  ont.raw.counts.per.cell.per.window = lapply(define.windows, function(x){
#    
#    meta.temp = m.meta[c(x[1]:x[2]),]
#    shared.cells <- intersect(rownames(meta.temp), rownames(ch.ont.samples.integrated@meta.data))
#    message(x[1],":", x[2])
#    
#    counts <- as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[gene,shared.cells])
#    colnames(counts)[1] <- c("gene.counts")
#    counts$total <- colSums(as.data.frame(ch.ont.samples.integrated@assays$RNA@counts[,shared.cells]))
#    counts$nCount_RNA <- ch.ont.samples.integrated@meta.data[shared.cells,]$nCount_RNA
#    counts$norm.counts <- (counts$gene.counts/counts$total) * 10000
#    counts$ln.norm.counts <- log1p(counts$norm.counts)
#    counts$seurat.norm.counts <- ch.ont.samples.integrated@assays$RNA@data[gene,shared.cells]
#    counts$genotype <- ch.ont.samples.integrated@meta.data[shared.cells,]$Genotype_1WT
#    
    #counts$used <- 0
    #counts[which(counts$nCount_RNA >= read.count.threshold),]$used <- 1
    
#    out = counts
#    return(out)
#    
#  })
#  
#  ont.raw.counts.per.cell.per.window.per.gene[[gene]] <- ont.raw.counts.per.cell.per.window
#}


##------------------------------- Save all output ------------------------------- ##

#save(pseudotime.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.celltype.unique.bucket.100.Rdata")
#save(pseudotime.mean, ct.fractions, ct.fractions.2, win.mutpsi, ont.mut.expression.per.window, ont.mut.raw.counts.per.window, ont.raw.counts.per.cell.per.window.per.gene, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.info.sig.in.any.dPSI2.bucket60.onlyMUT.Rdata")
save(pseudotime.mean, ct.fractions, ct.fractions.2, win.mutpsi,win.ct.counts, ont.mut.expression.per.window, ont.mut.raw.counts.per.window, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/CH.output/CH259.CH305.merged.info.sig.in.any.dPSI2.bucket60.onlyMUT.wCTcounts.Rdata")


