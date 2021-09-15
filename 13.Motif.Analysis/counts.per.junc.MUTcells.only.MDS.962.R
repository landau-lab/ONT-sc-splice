##------------------------------- load in needed libraries------------------------------- ##
library(tidyverse)
library(data.table)
library(Seurat)

##Motif Analysis
##Generating counts per Pseudotime window for MDS 

##Using the Object grab the HSPC, IMP, MEP, EP cells, rank by pseudotime and sperate into quantiles.

##----------------- Setting important variables   -----------------##

#seurat.object.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/mds.sample.integrated.p4.5.6.in.progress_2_withPseudotime.RData"

junction.meta.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/MDS962/output_files/leafcutter_outputs/MDS962_output/MDS962_all.introns.info.w.primary.annotations.txt"

junction.count.file.path <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/MDS962/output_files/leafcutter_outputs/MDS962_output/MDS962_perind_numbers.counts.txt"

genotyping.table.file.path <- "/gpfs/commons/home/pchamely/SequencingData/MDS_dataset_06_06_21/MDS962.DIG.genotyping.txt"

genotype.column <- "Genotype_Approx_Match_No_Gene_Thresh" #to take from genotyping table 

##----------------- Load in Seurat Object   -----------------##

#print("Loading in seurat object")
#load(seurat.object.path) #seurat
#mds.samples.integrated <- seurat

##----------------- Load in counts matrix -----------------##

print("Loading count matrix")
m <- fread(junction.count.file.path, data.table=FALSE)
rownames(m) <- m$intron_junction
m <- m[,-1]

##----------------- Grab the genotyped cells that we want to use (only MUT) -----------------##

genotyping.table <- read.table(genotyping.table.file.path)
#cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] =="MUT"),])
cells.to.use <- rownames(genotyping.table[which(genotyping.table[,genotype.column] %in% c("MUT","WT")),])

##----------------- Load the metadata & remove NA's from pseudotime column -----------------##

#cell.meta = as.data.frame(mds.samples.integrated@meta.data[cells.to.use,])
#cell.meta <- cell.meta[!is.na(cell.meta[,pseudotime.column]),]
shared.cells <- intersect(colnames(m), cells.to.use)

print("number of cell in long read data")
ncol(m)
print("number of mut genotyped cell in short read data")
length(cells.to.use)
print("number of shared cells")
length(shared.cells)

junction.meta <- read.csv(junction.meta.file.path)
junction.meta$intron_junction <- paste0(junction.meta$chr,":",junction.meta$start,":",junction.meta$end,":",junction.meta$strand)
junction.meta <- junction.meta[,-1]

print("num rows before removing duplicates")
nrow(junction.meta)
junction.meta <- junction.meta[!duplicated(junction.meta$intron_junction),]
print("num rows after removing duplicates")
nrow(junction.meta)

rownames(junction.meta) <- junction.meta$intron_junction
shared.junctions <- intersect(rownames(m), rownames(junction.meta))

#Subset the cell.metadata, junction count matrix and junction.metadata to have the same cells and junctions
genotyping.table.sub <- genotyping.table[shared.cells,]
m.sub = m[shared.junctions,shared.cells]
junction.meta.sub <- junction.meta[shared.junctions,]


##------------------------------- Strand adjust the metadata ------------------------------- ##

junction.meta.sub$startClass <- paste(junction.meta.sub$startClass)
junction.meta.sub$endClass <- paste(junction.meta.sub$endClass)


junction.meta.sub$fivep_class = 'NA'
junction.meta.sub$threep_class = 'NA'

junction.meta.sub = junction.meta.sub %>% mutate(fivep_class = ifelse(strand_1=="+", startClass, endClass),
                                                       threep_class = ifelse(strand_1 == "+", endClass, startClass))
print("Taking strandedness into account")

#Taking strandedness into account 
junction.meta.sub$fivep_distance <- NA
junction.meta.sub$threep_distance <- NA
junction.meta.sub[which(junction.meta.sub$strand_1 == "+"), ]$fivep_distance <- junction.meta.sub[which(junction.meta.sub$strand_1 == "+"), ]$startDistance
junction.meta.sub[which(junction.meta.sub$strand_1 == "+"), ]$threep_distance <- junction.meta.sub[which(junction.meta.sub$strand_1 == "+"), ]$endDistance
junction.meta.sub[which(junction.meta.sub$strand_1 == "-"),]$fivep_distance <- (junction.meta.sub[which(junction.meta.sub$strand_1 == "-"), ]$endDistance)*-1
junction.meta.sub[which(junction.meta.sub$strand_1 == "-"),]$threep_distance <- (junction.meta.sub[which(junction.meta.sub$strand_1 == "-"), ]$startDistance)*-1
print("Strand adjustment done")
print(dim(junction.meta.sub))
#subset metadata by those that have strand information and filter counts information for only those that have strand information 

junction.meta.sub$intron_junction = paste(junction.meta.sub$chr, junction.meta.sub$start, junction.meta.sub$end, junction.meta.sub$strand_1, sep = ":")

junction.meta.sub = junction.meta.sub %>% filter(strand_1 != 'NA')



junction.meta.sub$Final_Verdict = NA
junction.meta.sub[which(junction.meta.sub$startClass == "main" & junction.meta.sub$endClass =="main"),"Final_Verdict"] = "Canonical"
junction.meta.sub[which(((junction.meta.sub$startClass == "main" & junction.meta.sub$endClass =="not_main_3_prime") | (junction.meta.sub$startClass == "not_main_3_prime" & junction.meta.sub$endClass =="main")) &
                   (junction.meta.sub$threep_distance > (-100) & junction.meta.sub$threep_distance < 0)),"Final_Verdict"] = "Cryptic_threeprime"
junction.meta.sub[which(((junction.meta.sub$startClass == "main" & junction.meta.sub$endClass =="not_main_5_prime") | (junction.meta.sub$startClass == "not_main_5_prime" & junction.meta.sub$endClass =="main")) &
                   (junction.meta.sub$fivep_distance > (-100) & junction.meta.sub$fivep_distance < 0)),"Final_Verdict"] = "Cryptic_fiveprime"
junction.meta.sub[which(junction.meta.sub$startClass == "not_main_3_prime" & junction.meta.sub$endClass =="not_main_5_prime"),"Final_Verdict"] = "cryptic_unanchored"
junction.meta.sub[which(junction.meta.sub$startClass == "not_main_5_prime" & junction.meta.sub$endClass =="not_main_3_prime"),"Final_Verdict"] = "cryptic_unanchored"
junction.meta.sub[which(((junction.meta.sub$startClass == "not_main_3_prime" & junction.meta.sub$end == "main") | (junction.meta.sub$startClass == "main" & junction.meta.sub$endClass =="not_main_3_prime")) &
                   (junction.meta.sub$threep_distance < (-100) | junction.meta.sub$threep_distance > 0)), "Final_Verdict"] = "Alternative_threeprime"
junction.meta.sub[which(((junction.meta.sub$startClass == "not_main_5_prime" & junction.meta.sub$end == "main") | (junction.meta.sub$startClass == "main" & junction.meta.sub$endClass =="not_main_5_prime")) &
                   (junction.meta.sub$fivep_distance < (-100) | junction.meta.sub$fivep_distance > 0)), "Final_Verdict"] = "Alternative_fiveprime"

##------------------------------- Order cells in pseudotime and group by quantile ------------------------------- ##

#ord = rownames(cell.meta.sub)[order(cell.meta.sub[,pseudotime.column],decreasing = F)]
#cell.meta.sub = cell.meta.sub[ord,]
#m.sub = m.sub[,ord]

#num.quantiles <- 4
#cell.meta.sub$quantile <- ntile(cell.meta.sub[,pseudotime.column], num.quantiles)
#table(cell.meta.sub$quantile)


##------------------------------- Get mut.psi info for each quantile and add to table ------------------------------- ##

#for (i in 1:num.quantiles){
#  print(i)
#  #i <- 1
#  quantile.cells <- rownames(cell.meta.sub[which(cell.meta.sub$quantile == i),])
#  col.name <- paste0("Q",i,".junc.cov")
#  junction.meta.sub[,col.name] <- rowSums(m.sub[,quantile.cells])
#}


##------------------------------- Get mut cell counts for each junction ------------------------------- ##

#m.sub = m[shared.junctions,shared.cells]
#junction.meta.sub <- junction.meta[shared.junctions,]

#junction.meta.sub$mut.cells.junc.cov <- rowSums(m.sub)

mut.cells <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] == "MUT"),])
length(mut.cells)
junction.meta.sub$mut.cells.junc.cov <- rowSums(m.sub[,mut.cells])

wt.cells <- rownames(genotyping.table.sub[which(genotyping.table.sub[,genotype.column] == "WT"),])
length(wt.cells)
junction.meta.sub$wt.cells.junc.cov <- rowSums(m.sub[,wt.cells])


##------------------------------- Save results ------------------------------- ##

write.table(junction.meta.sub, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/MDS962.N626D.junction.info.w.mut.wt.cell.counts.full.txt", quote = T, sep=",",row.names = TRUE)

junction.meta.final <- junction.meta.sub[,c("intron_junction","clusterID_5p","gene","Final_Verdict" ,"fivep_distance","threep_distance","mut.cells.junc.cov", "wt.cells.junc.cov")]
colnames(junction.meta.final) <- c("intron.junction","clusterID","gene","final.junction.verdict","fivep.distance","threep.distance","mut.cells.junc.cov","wt.cells.junc.cov" )

write.table(junction.meta.final, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/13.Motif.Analysis/junction.matricies/MDS962.N626D.junction.info.w.mut.wt.cell.counts.txt", quote = F)

