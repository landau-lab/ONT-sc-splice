
##------------------------------- load in needed libraries------------------------------- ##

library(tidyverse)
library(data.table)

##------------------------------- File paths and arguments ------------------------------- ##

#junction.of.interest <- "chr1:67424977:67425082:-" 
#gene <- "SERBP1"

junction.of.interest <- "chr7:80672040:80672754:+" 
gene <- "CD36"

intron.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.intron.metadata.txt"
junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.count.matrix.txt"

##------------------------------- Identify cluster of junctions from intron metadata file ------------------------------- ##

print("Grapping clusterID and junction set")
intron.meta <- read.table(intron.meta.file.path, sep="\t", header = T)
intron.meta$intron.junction <- paste0(intron.meta$chr,":",intron.meta$start,":",intron.meta$end,":",intron.meta$strand)
clusterID <- paste(intron.meta[which(intron.meta$intron.junction == junction.of.interest),]$clusterID_5p)
intron.meta.sub <- intron.meta[which(intron.meta$clusterID_5p == clusterID),]
rownames(intron.meta.sub) <- intron.meta.sub$intron.junction


##------------------------------- Subset count matrix for cluster of junctions ------------------------------- ##

print("Loading count matrix")
count.matrix <- fread(junction.count.file.path, data.table=FALSE)
rownames(count.matrix) <- count.matrix$intron_junction
count.matrix <- count.matrix[-1]

shared.junctions <- intersect(rownames(intron.meta.sub), rownames(count.matrix))
print("Number of junctions in cluster")
length(rownames(intron.meta.sub))

print("Number of shared junctions")
length(shared.junctions)

count.matrix.sub <- count.matrix[shared.junctions,]


##------------------------------- Save new sub count matrix output ------------------------------- ##


write.table(count.matrix.sub, file=paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/16.Pseudotime.Window/IterativeDownsampling/cluster.count.matrix.subsets/",gene,".",clusterID,".count.matrix.sub.txt"))

