##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##Motif Analysis
##Generating counts per Pseudotime window for MDS 

##Using the Object grab the HSPC, IMP, MEP, EP cells, rank by pseudotime and sperate into quantiles.

##----------------- Setting important variables   -----------------##

seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ch.integrated.Rdata"

#junction.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/1.ind_patient_logOR/100k_perm/CH259/logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt"

junction.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/CH259.K666N.alt3p.junction.info.w.quantile.and.ct.counts.mutcells.only.txt"

junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_perind_numbers.counts.txt"

genotyping.table.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt"

#pseudotime.column <- "pseudotime_monocle3" #to take from seurat object
genotype.column <- "Genotype_2UMI" #to take from genotyping table 
cluster.column = "Cell.Assignment_2" #to take from genotyping table
#cluster.types = c("HSPC","IMP","MEP","EP")

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
colnames(m) <- paste0(colnames(m),"_1")
colnames(m)[1:5]


##----------------- Grab the genotyped cells that we want to use (only MUT) -----------------##

genotyping.table <- read.table(genotyping.table.file.path)
#cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] %in% c("MUT","WT") & !is.na(genotyping.table[,cluster.column]) ),])
cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] == "MUT" & !is.na(genotyping.table[,cluster.column]) ),])
##----------------- Load the metadata & remove NA's from pseudotime column -----------------##

shared.cells <- intersect(cells.to.use, rownames(ch.samples.integrated@meta.data))
length(shared.cells)
length(cells.to.use)
cell.meta = as.data.frame(ch.samples.integrated@meta.data[shared.cells,])
cell.meta[1:5,1:5]

#cell.meta <- cell.meta[!is.na(cell.meta[,pseudotime.column]),]
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
junction.meta <- junction.meta[!duplicated(junction.meta$intron.junction),]
print("num rows after removing duplicates")
nrow(junction.meta)

print("now adding junctions as row.names")
rownames(junction.meta) <- junction.meta$intron.junction
print("finished adding junctions as row.names")
shared.junctions <- intersect(rownames(m), rownames(junction.meta))
print("number of shared junctions")
length(shared.junctions)


#Subset the cell.metadata, junction count matrix and junction.metadata to have the same cells and junctions
cell.meta.sub = cell.meta[shared.cells,]
m.sub = m[shared.junctions,shared.cells]
junction.meta.sub <- junction.meta[shared.junctions,]

##------------------------------- Get mut cell counts for each junction ------------------------------- ##

#m.sub = m[shared.junctions,shared.cells]
#junction.meta.sub <- junction.meta[shared.junctions,]

mut.cells <- rownames(genotyping.table[which(genotyping.table[,genotype.column] == "MUT"),])
length(mut.cells)
shared.mut <- intersect(mut.cells,colnames(m.sub))
print("num.shared.mut.cells")
length(shared.mut)
junction.meta.sub$bulk.mut.cells.junc.cov <- rowSums(m.sub[,shared.mut])

#wt.cells <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] == "WT"),])
#length(wt.cells)
#shared.wt <- intersect(wt.cells,colnames(m.sub))
#print("num.shared.wt.cells")
#length(shared.wt)
#junction.meta.sub$wt.cells.junc.cov <- rowSums(m.sub[,shared.wt])



##------------------------------- Save results ------------------------------- ##

write.table(junction.meta.sub, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/CH259.K666N.alt3p.junction.info.w.quantile.and.ct.and.bulk.counts.mutcells.only.txt", quote = F)


