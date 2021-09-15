
##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##------------------------------- Generate a list of cryptic junctions to look into ------------------------------- ##

load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/junctions.to.use.MDS.noSigThresh.ct.sub.Rdata")

print("junctions")
junctions.of.interest
print("junctions with gene names")
junctions.of.interest.final
print("gene list")
gene.list

##------------------------------- File paths and arguments ------------------------------- ##
seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/mds.sample.integrated.p4.5.6.in.progress_2_withPseudotime.RData"
intron.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.intron.metadata.txt"
junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.count.matrix.txt"
genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.genotype.info.with.cell.type.txt"


all.cite.markers <- c("ADT-CD71","ADT-CD36","ADT-CD99","ADT-CD38")
main.cite.marker <- "ADT-CD71"
#pseudotime.column <- "pseudotime_monocole"
genotype.column <- "Genotype_2UMI_mutUMIFrac0.2" #to take from genotyping table 
step  <- 300
cluster.column = "Cell.Assignment_adj" 
#cluster.types = c("HSPC","IMP","NP","MkP","MEP","EP", "Mat_Mono")
cluster.types = c("HSPC","MkP","MEP","EP")

##------------------------------- load in seurat object, intron metadata junction count matrix and genotyping.info files ------------------------------- ##

load(seurat.object.path) #seurat
mds.samples.integrated <- seurat

intron.meta <- read.table(intron.meta.file.path, sep="\t", header = T)
intron.meta$intron.junction <- paste0(intron.meta$chr,":",intron.meta$start,":",intron.meta$end,":",intron.meta$strand)
clusterIDs <- paste(intron.meta[which(intron.meta$intron.junction %in% junctions.of.interest),]$clusterID_5p)
intron.meta.sub <- intron.meta[which(intron.meta$clusterID_5p %in% clusterIDs),]
rownames(intron.meta.sub) <- intron.meta.sub$intron.junction

print("Loading count matrix")
#m <- read.table(junction.count.file.path, header = T)
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[,-1]

print("Done loading count matrix")
m[1:5,1:5]

genotyping.table <- read.table(genotyping.table.file.path)


##------------------------------- Add in the Cite expression to metadata and generate metadata table (subsetting for only genotyped cells) & remove NA's from main.cite.marker column ------------------------------- ##

DefaultAssay(mds.samples.integrated) <- "ADT"
ADT_expression <- FetchData(mds.samples.integrated, vars= all.cite.markers)
mds.samples.integrated <- AddMetaData(mds.samples.integrated, ADT_expression)
all.cite.marker.columns <- gsub("-", ".", all.cite.markers)
main.cite.marker.column <- gsub("-", ".", main.cite.marker)
DefaultAssay(mds.samples.integrated) <- "RNA"


meta = as.data.frame(mds.samples.integrated@meta.data[which(mds.samples.integrated@meta.data[,cluster.column] %in% cluster.types & mds.samples.integrated@meta.data$Genotype_2UMI == "MUT"),])

meta <- meta[!is.na(meta[,main.cite.marker.column]),]

#Subset the metadata and junction count matrix to have the same cells (and remove any NAs)
shared.cells <- intersect(colnames(m), rownames(meta))
print("number of cells in long read data")
ncol(m)
print("number of cells in short read data")
nrow(meta)
print("number of shared cells")
length(shared.cells)

m.sub = m[,shared.cells]
m.meta = meta[shared.cells,]

print("number of cells from each patient")
table(m.meta$orig.ident)

ncol(m.sub)
nrow(m.meta)

##------------------------------- Order cells in by cite marker expression ------------------------------- ##


ord = rownames(m.meta)[order(m.meta[,main.cite.marker.column],decreasing = F)]
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


##------------------------------- Get mean expression value and cell type fractions per window ------------------------------- ##

cite.expression.mean = lapply(define.windows, function(x){
  meta.temp = m.meta[c(x[1]:x[2]),]
  colMeans(meta.temp[,all.cite.marker.columns])
})


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
    
    ct.mut.cells <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] =="MUT" & genotyping.table.sub[,"Cell.Assignment"] == ct),])
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
        
       if (ct =="MkP"){
          message("Adding MkP counts")
          intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(MkP.cluster.cov.mut=sum(MkP.mut.bulk.counts))
          intron.meta.sub.temp <- intron.meta.sub.temp %>% group_by(clusterID_5p) %>% mutate(MkP.mut.psi = (MkP.mut.bulk.counts/MkP.cluster.cov.mut)*100)
          rownames(intron.meta.sub.temp) <- rownames(intron.meta.sub)
          #print(intron.meta.sub.temp[130:150,c("MkP.mut.bulk.counts","MkP.cluster.cov.mut","MkP.mut.psi")])
          
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

load("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ONT.seruat.objects/mds.ont.samples.integrated.rna.normscale.adjusted.rownames.Rdata")
ont.samples.integrated <- mds.ont.samples.integrated

length(gene.list)
gene.list <- unique(intersect(rownames(ont.samples.integrated@assays$RNA@data), gene.list))
length(gene.list)

if (length(setdiff(gene.list, rownames(ont.samples.integrated@assays$RNA@data))) >0 ){
  message(("gene's not in Gene expression matrix - probably due to differences in annotation"))
  message(setdiff(gene.list, rownames(ont.samples.integrated@assays$RNA@data)))
}

ont.mut.expression.per.window = lapply(define.windows, function(x){
  
  meta.temp = m.meta[c(x[1]:x[2]),]
  shared.cells <- intersect(rownames(meta.temp), ont.samples.integrated@meta.data$adjusted.rownames)
  message(length(shared.cells))
  mut <- intersect(shared.cells, rownames(genotyping.table[which(genotyping.table[,genotype.column] =="MUT"),]))
  cell.to.use <- rownames(ont.samples.integrated@meta.data[which(ont.samples.integrated@meta.data$adjusted.rownames %in% mut),])
  message(length(cell.to.use))
  
  message(paste(length(mut),"MUT cells"))
  message(x[1],":", x[2])
  
  ont.samples.integrated@meta.data$window <- 0
  ont.samples.integrated@meta.data[cell.to.use,]$window <- 1
  
  Idents(ont.samples.integrated) <- "window"
  expression.values <- AverageExpression(ont.samples.integrated, assay = "RNA" ,features = gene.list, slot="data")
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
  shared.cells <- intersect(rownames(meta.temp), ont.samples.integrated@meta.data$adjusted.rownames)
  mut <- intersect(shared.cells, rownames(genotyping.table[which(genotyping.table[,genotype.column] =="MUT"),]))
  cell.to.use <- rownames(ont.samples.integrated@meta.data[which(ont.samples.integrated@meta.data$adjusted.rownames %in% mut),])
  length(cell.to.use)
  
  message(paste(length(mut),"MUT cells"))
  message(x[1],":", x[2])
  bulk.counts <- as.data.frame(rowSums(as.data.frame(ont.samples.integrated@assays$RNA@counts[gene.list,cell.to.use])))
  colnames(bulk.counts) <- c("counts")
  out = bulk.counts
  return(out)
  
})

ont.mut.raw.counts.per.window <- as.data.frame(ont.mut.raw.counts.per.window)
colnames(ont.mut.raw.counts.per.window) <- as.character(1:ncol(ont.mut.raw.counts.per.window))


#ont.raw.counts.per.cell.per.window.per.gene <- list()

#for (gene in gene.list){
  
#  ont.raw.counts.per.cell.per.window = lapply(define.windows, function(x){
#    
#    meta.temp = m.meta[c(x[1]:x[2]),]
#    shared.cells <- intersect(rownames(meta.temp), ont.samples.integrated@meta.data$adjusted.rownames)
#    message(x[1],":", x[2])
#    cell.to.use <- rownames(ont.samples.integrated@meta.data[which(ont.samples.integrated@meta.data$adjusted.rownames %in% shared.cells),])
#    
#    counts <- as.data.frame(ont.samples.integrated@assays$RNA@counts[gene,cell.to.use])
#    colnames(counts)[1] <- c("gene.counts")
#    counts$total <- colSums(as.data.frame(ont.samples.integrated@assays$RNA@counts[,cell.to.use]))
#    counts$nCount_RNA <- ont.samples.integrated@meta.data[cell.to.use,]$nCount_RNA
#    counts$norm.counts <- (counts$gene.counts/counts$total) * 10000
#    counts$ln.norm.counts <- log1p(counts$norm.counts)
#    counts$seurat.norm.counts <- ont.samples.integrated@assays$RNA@data[gene,cell.to.use]
#    counts$genotype <- genotyping.table[shared.cells,genotype.column]
#    
#    out = counts
#    return(out)
#    
#  })
#  
#  ont.raw.counts.per.cell.per.window.per.gene[[gene]] <- ont.raw.counts.per.cell.per.window
#}


##------------------------------- Save all output ------------------------------- ##

#save(cite.expression.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.mut.and.wt.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window,ont.raw.counts.per.cell.per.window.per.gene, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/MDS.output/MDSP5.MDSP6.merged.info.sig.in.any.bucket.300.usingCITE.ct.all.Rdata")
#save(cite.expression.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.mut.and.wt.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window,ont.raw.counts.per.cell.per.window.per.gene, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/MDS.output/MDSP5.MDSP6.merged.info.sig.in.any.bucket.400.usingCITE.ct.all.Rdata")


#save(cite.expression.mean, ct.fractions, ct.fractions.2, mut.fraction, win.dPSI, ont.expression.per.window, ont.mut.expression.per.window, ont.wt.expression.per.window, ont.mut.and.wt.raw.counts.per.window , ont.mut.raw.counts.per.window, ont.wt.raw.counts.per.window,ont.raw.counts.per.cell.per.window.per.gene, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/MDS.output/MDSP5.MDSP6.merged.info.sig.in.any.bucket.300.usingCITE.ct.sub.Rdata")

save(cite.expression.mean, ct.fractions, ct.fractions.2, win.mutpsi,win.ct.counts, ont.mut.expression.per.window, ont.mut.raw.counts.per.window, version=2, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/MDS.output/MDSP5.MDSP6.merged.info.no.sig.thresh.bucket.300.usingCITE.all.markers.ct.sub.onlyMut.wCTcounts.Rdata")


