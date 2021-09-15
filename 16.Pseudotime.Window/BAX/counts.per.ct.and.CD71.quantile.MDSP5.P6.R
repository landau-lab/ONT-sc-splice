##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##Motif Analysis
##Generating counts per Pseudotime window for MDS 

##Using the Object grab the HSPC, IMP, MEP, EP cells, rank by pseudotime and sperate into quantiles.

##----------------- Setting important variables   -----------------##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/mds.sample.integrated.p4.5.6.in.progress_2_withPseudotime.RData"

junction.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/MDS.untreated.combined.alt3p.junction.info.txt"

junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.count.matrix.txt"

#genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.genotype.info.with.cell.type.txt"
genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.genotype.info.with.cell.type.txt"

cite.marker <- "ADT-CD71"
genotype.column <- "Genotype_1UMI_mutUMIFrac0.2" 
cluster.column = "Cell.Assignment"
min.num.umis <- 1
cluster.types = c("HSPC","IMP","MEP","EP")
patient.list = c("MDS_P6")

##----------------- Load in Seurat Object   -----------------##

print("Loading in seurat object")
load(seurat.object.path) #seurat
mds.samples.integrated <- seurat

##----------------- Load in counts matrix -----------------##

print("Loading count matrix")
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[,-1]
colnames(m)[1:5]

##----------------- Grab the genotyped cells that we want to use (only MUT) -----------------##

genotyping.table <- read.table(genotyping.table.file.path)
cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] %in% c("MUT","WT") & genotyping.table[,cluster.column] %in% cluster.types & genotyping.table[,"total_10x"] >= min.num.umis),])
mut.cells <- rownames(genotyping.table[which(genotyping.table[,genotype.column] == "MUT" & genotyping.table[,cluster.column] %in% cluster.types & genotyping.table[,"total_10x"] >= min.num.umis),])
wt.cells <- rownames(genotyping.table[which(genotyping.table[,genotype.column] == "WT" & genotyping.table[,cluster.column] %in% cluster.types & genotyping.table[,"total_10x"] >= min.num.umis),])
length(cells.to.use)
length(mut.cells)
length(wt.cells)

##----------------- Load the metadata, add CITE expressoin & remove NA's from pseudotime column -----------------##

DefaultAssay(mds.samples.integrated) <- "ADT"
ADT_expression <- FetchData(mds.samples.integrated, vars= cite.marker)
mds.samples.integrated <- AddMetaData(mds.samples.integrated, ADT_expression)
cite.marker.column <- gsub("-", ".", cite.marker)
DefaultAssay(mds.samples.integrated) <- "RNA"


cell.meta = as.data.frame(mds.samples.integrated@meta.data[cells.to.use,])
print("checking patient ID numbers before subsetting")
table(cell.meta$orig.ident)
cell.meta = cell.meta[which(cell.meta$orig.ident %in% patient.list),]
print("checking patient ID numbers after subsetting")
table(cell.meta$orig.ident)
cell.meta[1:5,1:5]

cell.meta <- cell.meta[!is.na(cell.meta[,cite.marker.column]),]
shared.cells <- intersect(colnames(m), rownames(cell.meta))

print("number of cell in long read data")
ncol(m)
print("number of cell in short read data")
nrow(cell.meta)
print("number of shared cells")
length(shared.cells)

junction.meta <- read.table(junction.meta.file.path)
#junction.meta$intron.junction <- paste0(junction.meta$chr,":",junction.meta$start,":",junction.meta$end,":",junction.meta$strand)
#junction.meta <- junction.meta[,-1]
print("num rows before removing duplicates")
nrow(junction.meta)
junction.meta <- junction.meta[!duplicated(junction.meta$intron.junction),]
print("num rows after removing duplicates")
nrow(junction.meta)

print("now adding junctions as row.names")
rownames(junction.meta) <- junction.meta$intron.junction
rownames(m)[1:5]
rownames(junction.meta)[1:5]
print("finished adding junctions as row.names")
shared.junctions <- intersect(rownames(m), rownames(junction.meta))
print("number of shared junctions")
length(shared.junctions)

#Subset the cell.metadata, junction count matrix and junction.metadata to have the same cells and junctions
cell.meta.sub = cell.meta[shared.cells,]
m.sub = m[shared.junctions,shared.cells]
junction.meta.sub <- junction.meta[shared.junctions,]

print("checking patient ID numbers again")
table(cell.meta.sub$orig.ident)



##------------------------------- Order cells in pseudotime and group by quantile ------------------------------- ##

ord = rownames(cell.meta.sub)[order(cell.meta.sub[,cite.marker.column],decreasing = F)]
cell.meta.sub = cell.meta.sub[ord,]
m.sub = m.sub[,ord]

num.quantiles <- 4
cell.meta.sub$quantile <- ntile(cell.meta.sub[,cite.marker.column], num.quantiles)
table(cell.meta.sub$quantile)

##------------------------------- Get mut.psi and wt.psi overall and add to table ------------------------------- ##

junction.meta.sub[,"bulk.mut.junc.cov"] <- rowSums(m.sub[,intersect(rownames(cell.meta.sub) ,mut.cells)]) 
junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(bulk.mut.cluster.cov = sum(bulk.mut.junc.cov))
junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(bulk.mut.bulk.psi = (bulk.mut.junc.cov/bulk.mut.cluster.cov)*100)


junction.meta.sub[,"bulk.wt.junc.cov"] <- rowSums(m.sub[,intersect(rownames(cell.meta.sub) ,wt.cells)])
junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(bulk.wt.cluster.cov = sum(bulk.wt.junc.cov))
junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(bulk.wt.bulk.psi = (bulk.wt.junc.cov/bulk.wt.cluster.cov)*100)



##------------------------------- Get mut.psi info for each quantile and add to table ------------------------------- ##

#junction.meta.sub[which(junction.meta.sub$gene =="BAX"),]

qunatile.count.info <- as.data.frame(matrix(ncol=2,nrow=num.quantiles))
colnames(qunatile.count.info) <- c("num.mut.cells","num.wt.cells")
rownames(qunatile.count.info) <- 1:num.quantiles

for (i in 1:num.quantiles){
  #print(i)
  quantile.mut.cells <- intersect(rownames(cell.meta.sub[which(cell.meta.sub$quantile == i ),]),mut.cells)
  qunatile.count.info[i, "num.mut.cells"] <- length(quantile.mut.cells)
  col.name <- paste0("Q",i,".mut.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,quantile.mut.cells]) 
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0("Q",i,".mut.cluster.cov") := sum(!!as.name(paste0(col.name))))
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0("Q",i,".mut.bulk.psi") := (!!as.name(paste0(col.name))/!!as.name(paste0("Q",i,".mut.cluster.cov")))*100)
  
  #junction.meta.sub[which(junction.meta.sub$gene =="BAX"),]
  
  quantile.wt.cells <- intersect(rownames(cell.meta.sub[which(cell.meta.sub$quantile == i ),]),wt.cells)
  qunatile.count.info[i, "num.wt.cells"] <- length(quantile.wt.cells)
  col.name <- paste0("Q",i,".wt.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,quantile.wt.cells])
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0("Q",i,".wt.cluster.cov") := sum(!!as.name(paste0(col.name))))
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0("Q",i,".wt.bulk.psi") := (!!as.name(paste0(col.name))/!!as.name(paste0("Q",i,".wt.cluster.cov")))*100)
  
  
}


##------------------------------- Get mut.psi info for each cell.type and add to table ------------------------------- ##

cell.count.info <- as.data.frame(matrix(ncol=2,nrow=length(cluster.types)))
colnames(cell.count.info) <- c("num.mut.cells","num.wt.cells")
rownames(cell.count.info) <- cluster.types

for (ct in cluster.types){
  print(ct)
  ct.mut.cells <- intersect(rownames(cell.meta.sub[which(cell.meta.sub$Cell.Assignment_adj == ct ),]),mut.cells)
  cell.count.info[ct, "num.mut.cells"] <- length(ct.mut.cells)
  col.name <- paste0(ct,".mut.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,ct.mut.cells]) 
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0(ct,".mut.cluster.cov") := sum(!!as.name(paste0(col.name))))
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0(ct,".mut.bulk.psi") := (!!as.name(paste0(col.name))/!!as.name(paste0(ct,".mut.cluster.cov")))*100)
  
  ct.wt.cells <- intersect(rownames(cell.meta.sub[which(cell.meta.sub$Cell.Assignment_adj == ct ),]),wt.cells)
  cell.count.info[ct, "num.wt.cells"] <- length(ct.wt.cells)
  col.name <- paste0(ct,".wt.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,ct.wt.cells])
  
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0(ct,".wt.cluster.cov") := sum(!!as.name(paste0(col.name))))
  junction.meta.sub <- junction.meta.sub %>% group_by(gene, clusterID) %>% mutate(!!paste0(ct,".wt.bulk.psi") := (!!as.name(paste0(col.name))/!!as.name(paste0(ct,".wt.cluster.cov")))*100)
  

}


print("num of mut and wt cells")

print("per quantile")
print(qunatile.count.info)


print("per.ct")
print(cell.count.info)



##------------------------------- Save results ------------------------------- ##

write.table(junction.meta.sub, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/BAX/MDS.alt3p.junction.info.w.quantile.and.ct.counts.minUMI.1.MDSP6.txt", quote = T,row.names = TRUE, na = "NA",sep=",")

