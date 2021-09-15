#!/usr/bin/env Rscript
#options(echo=TRUE)
### ONT - Generating merged bulk PSI scores

##--------------------------- Load in necessary libraries ---------------------------##

library(tibble)
library(tidyverse)
library(optparse)

##--------------------------- Allow for the input files to be arguments and save them as variables ---------------------------##

option_parser=OptionParser(
  usage="%prog [options] sample.list intron.metadata.with.bulkPSI.1.tsv intron.metadata.with.bulkPSI.2.tsv output.file"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 4)

sample.list <- unlist(strsplit(parsed_args$args[1], ",")[[1]])

print(sample.list[1])
print(sample.list[2])

intron.meta.data.files <- list()
intron.meta.data.files[[sample.list[1]]] <- parsed_args$args[2]
intron.meta.data.files[[sample.list[2]]] <- parsed_args$args[3]

print(intron.meta.data.files[[sample.list[1]]])
print(intron.meta.data.files[[sample.list[2]]])

output.file <- parsed_args$args[4]


##--------------------------- Merge Patient counts to get combined Bulk PSI - using the *intron.metadata.with.bulkPSI.ONT.tsv files ---------------------------##

meta.data.table.list <- list()

##Read in intron.metadat.table1
for (sample in sample.list){
  meta.data.table.list[[sample]] <- read.table(intron.meta.data.files[[sample]], sep="\t")
  rownames(meta.data.table.list[[sample]]) <- meta.data.table.list[[sample]]$intron.junction
}


##Find intersecting junctions and create new merged table
shared.junctions <- intersect(rownames(meta.data.table.list[[1]]) ,  rownames(meta.data.table.list[[2]]))
intron.meta.data.merged <- meta.data.table.list[[1]][shared.junctions,c("intron.junction","chr","start","end","strand","five_prime","three_prime","gene")]

##Merge the final counts into new table (both the individual junction coverage and the cluster coverage) 
for (sample in sample.list){
  intron.meta.data.merged[,paste0("total.junction.coverage.",sample)] <- meta.data.table.list[[sample]][shared.junctions,]$total.junction.coverage
  intron.meta.data.merged[,paste0("cluster.coverage.",sample)] <- meta.data.table.list[[sample]][shared.junctions,]$cluster.coverage
  intron.meta.data.merged[,paste0("bulk.psi.",sample)] <- meta.data.table.list[[sample]][shared.junctions,]$bulk.psi
}


intron.meta.data.merged[,"total.junction.coverage"] <- intron.meta.data.merged[,paste0("total.junction.coverage.",sample.list[1])] +  intron.meta.data.merged[,paste0("total.junction.coverage.",sample.list[2])]

intron.meta.data.merged[,"cluster.coverage"] <- intron.meta.data.merged[,paste0("cluster.coverage.",sample.list[1])] +  intron.meta.data.merged[,paste0("cluster.coverage.",sample.list[2])]

intron.meta.data.merged[,"bulk.psi"] <- (intron.meta.data.merged[,"total.junction.coverage"]/intron.meta.data.merged[,"cluster.coverage"])*100


##--------------------------- Save output table with merged bulk PSI scores ---------------------------##

write.table(intron.meta.data.merged, file=output.file, quote = F, sep="\t")


