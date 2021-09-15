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

genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.genotype.info.with.cell.type.txt"

pseudotime.column <- "pseudotime_monocole"
genotype.column <- "Genotype_2UMI_mutUMIFrac0.2" #to take from genotyping table 
cluster.column = "Cell.Assignment_adj" 
cluster.types = c("HSPC","IMP","MEP","EP")

##----------------- Load in Seurat Object   -----------------##

print("Loading in seurat object")
load(seurat.object.path) #seurat
mds.samples.integrated <- seurat

##----------------- Load in counts matrix -----------------##

print("Loading count matrix")
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[,-1]


##----------------- Grab the genotyped cells that we want to use (only MUT) -----------------##

genotyping.table <- read.table(genotyping.table.file.path)
cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] =="MUT" & genotyping.table$Cell.Assignment %in% cluster.types),])

##----------------- Load the metadata & remove NA's from pseudotime column -----------------##

cell.meta = as.data.frame(mds.samples.integrated@meta.data[cells.to.use,])
cell.meta <- cell.meta[!is.na(cell.meta[,pseudotime.column]),]
shared.cells <- intersect(colnames(m), rownames(cell.meta))

print("number of cell in long read data")
ncol(m)
print("number of cell in short read data")
nrow(cell.meta)
print("number of shared cells")
length(shared.cells)

junction.meta <- read.table(junction.meta.file.path)
print("num rows before removing duplicates")
nrow(junction.meta)
junction.meta <- junction.meta[!duplicated(junction.meta$intron.junction),]
print("num rows after removing duplicates")
nrow(junction.meta)

rownames(junction.meta) <- junction.meta$intron.junction
shared.junctions <- intersect(rownames(m), rownames(junction.meta))

#Subset the cell.metadata, junction count matrix and junction.metadata to have the same cells and junctions
cell.meta.sub = cell.meta[shared.cells,]
m.sub = m[shared.junctions,shared.cells]
junction.meta.sub <- junction.meta[shared.junctions,]


##------------------------------- Order cells in pseudotime and group by quantile ------------------------------- ##

ord = rownames(cell.meta.sub)[order(cell.meta.sub[,pseudotime.column],decreasing = F)]
cell.meta.sub = cell.meta.sub[ord,]
m.sub = m.sub[,ord]

num.quantiles <- 4
cell.meta.sub$quantile <- ntile(cell.meta.sub[,pseudotime.column], num.quantiles)
table(cell.meta.sub$quantile)


##------------------------------- Get mut.psi info for each quantile and add to table ------------------------------- ##

for (i in 1:num.quantiles){
  print(i)
  #i <- 1
  quantile.cells <- rownames(cell.meta.sub[which(cell.meta.sub$quantile == i),])
  col.name <- paste0("Q",i,".junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,quantile.cells])
}



##------------------------------- Save results ------------------------------- ##

write.table(junction.meta.sub, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/MDS.untreated.combined.alt3p.junction.info.w.quantile.counts.txt", quote = F)


