##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##Motif Analysis
##Generating counts per Pseudotime window for MDS 

##Using the Object grab the HSPC, IMP, MEP, EP cells, rank by pseudotime and sperate into quantiles.

##----------------- Setting important variables   -----------------##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ch.integrated.Rdata"

junction.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/1.ind_patient_logOR/100k_perm/CH305/logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt"

junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH305/output_files/leafcutter_outputs/CH305_output/CH305_perind_numbers.counts.txt"

genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_305.genotype.info.with.cell.type.txt"

pseudotime.column <- "pseudotime_monocle3" #to take from seurat object
genotype.column <- "Genotype_2UMI" #to take from genotyping table 
cluster.column = "Cell.Assignment_2" #to take from genotyping table
cluster.types = c("HSPC","IMP","MEP","EP")

##----------------- Load in Seurat Object   -----------------##

print("Loading in seurat object")
load(seurat.object.path) #samples.integrated
ch.samples.integrated <- samples.integrated

##----------------- Load in counts matrix -----------------##

print("Loading count matrix")
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[,-1]
colnames(m)[1:5]
colnames(m) <- paste0(colnames(m),"_2")
colnames(m)[1:5]


##----------------- Grab the genotyped cells that we want to use (only MUT) -----------------##

genotyping.table <- read.table(genotyping.table.file.path)
cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] =="MUT" & genotyping.table[,cluster.column] %in% cluster.types),])

##----------------- Load the metadata & remove NA's from pseudotime column -----------------##

cell.meta = as.data.frame(ch.samples.integrated@meta.data[cells.to.use,])
cell.meta[1:5,1:5]

cell.meta <- cell.meta[!is.na(cell.meta[,pseudotime.column]),]
shared.cells <- intersect(colnames(m), rownames(cell.meta))

print("number of cell in long read data")
ncol(m)
print("number of cell in short read data")
nrow(cell.meta)
print("number of shared cells")
length(shared.cells)

junction.meta <- read.table(junction.meta.file.path)
#junction.meta <- read.csv(junction.meta.file.path)
#junction.meta$intron.junction <- paste0(junction.meta$chr,":",junction.meta$start,":",junction.meta$end,":",junction.meta$strand)
#junction.meta <- junction.meta[,-1]
print("num rows before removing duplicates")
nrow(junction.meta)
junction.meta <- junction.meta[!duplicated(junction.meta$intron_junction),]
print("num rows after removing duplicates")
nrow(junction.meta)

print("now adding junctions as row.names")
rownames(junction.meta) <- junction.meta$intron_junction
print("finished adding junctions as row.names")
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
  col.name <- paste0("Q",i,".mut.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,quantile.cells])
}

##------------------------------- Get mut.psi info for each cell.type and add to table ------------------------------- ##

for (ct in cluster.types){
  print(ct)
  ct.cells <- rownames(cell.meta.sub[which(cell.meta.sub$CellType == ct),])
  col.name <- paste0(ct,".mut.junc.cov")
  junction.meta.sub[,col.name] <- rowSums(m.sub[,ct.cells])

}


##------------------------------- Save results ------------------------------- ##

write.table(junction.meta.sub, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/CH305.K700E.alt3p.junction.info.w.quantile.and.ct.counts.mutcells.only.full.txt", quote = T,row.names = TRUE, na = "NA",sep=",")

junction.meta.sub.final <- junction.meta.sub[,c("pvalue","intron_junction","five_prime_ID","gene","Final_Verdict","fivep_distance","threep_distance", "Q1.mut.junc.cov","Q2.mut.junc.cov","Q3.mut.junc.cov","Q4.mut.junc.cov","HSPC.mut.junc.cov","IMP.mut.junc.cov","MEP.mut.junc.cov","EP.mut.junc.cov")]


colnames(junction.meta.sub.final) <- c("pvalue","intron.junction","clusterID","gene","final.junction.verdict","fivep.distance","threep.distance","Q1.mut.junc.cov","Q2.mut.junc.cov","Q3.mut.junc.cov","Q4.mut.junc.cov","HSPC.mut.junc.cov","IMP.mut.junc.cov","MEP.mut.junc.cov","EP.mut.junc.cov")

write.table(junction.meta.sub.final, file=paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/CH305.K700E.alt3p.junction.info.w.quantile.and.ct.counts.mutcells.only.txt"), quote = F)
