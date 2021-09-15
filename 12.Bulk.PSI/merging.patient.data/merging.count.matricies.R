#!/usr/bin/env Rscript
#options(echo=TRUE)
### Merge the count matricies  

##--------------------------- Load in necessary libraries ---------------------------##

library(tidyverse)
library(optparse)
library(data.table)

##--------------------------- Allow for the input files to be arguments and save them as variables ---------------------------##

option_parser=OptionParser(
  usage="%prog [options] counts.matrix.list.file output.file.prefix"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 2)

count.matricies.paths <- parsed_args$args[1]
output.file.prefix <- parsed_args$args[2]

##--------------------------- Load in the .txt file with path to count matricies and convert to a list ---------------------------##

count.matricies <- scan(count.matricies.paths, what="", sep="\n")
print("number of count matricies being merged")
length(count.matricies)
print("count matricies being merged")
print(count.matricies)

#count.matrix1 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH259/output_files/leafcutter_outputs/CH259_output/CH259_perind_numbers.counts.txt", header = T)
#count.matrix2 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH305/output_files/leafcutter_outputs/CH305_output/CH305_perind_numbers.counts.txt", header = T)

count.matricies.list <- list()
for (i in 1:length(count.matricies)){
  print(i)
  count.matricies.list[[i]] <- fread(count.matricies[i], data.table=FALSE)
  ncols <- ncol(count.matricies.list[[i]])
  var <- (i+1)
  colnames(count.matricies.list[[i]])[2:ncols] <- paste0(colnames(count.matricies.list[[i]])[2:ncols],"_",var)
}

print("done reading in count maricies")
count.matricies.list[[1]][1:10,1:10]
count.matricies.list[[2]][1:10,1:10]
count.matricies.list[[3]][1:10,1:10]

##--------------------------- Merge the count matricies and convert "NAs" to 0's ---------------------------##

count.matricies.merged <- Reduce(function(d1, d2) merge(d1, d2, all.x = T, all.y = T), count.matricies.list)
count.matricies.merged[is.na(count.matricies.merged)] <- 0

print("number of cell for each")
print(ncol(count.matricies.list[[1]]))
print(ncol(count.matricies.list[[2]]))
print("number of cells merged")
print(ncol(count.matricies.merged))

print("number of rows for each")
print(nrow(count.matricies.list[[1]]))
print(nrow(count.matricies.list[[2]]))
print("number of rows with intersect")
print(length(union(count.matricies.list[[1]]$intron_junction, count.matricies.list[[2]]$intron_junction)))
print("number of rows merged")
print(nrow(count.matricies.merged))

##--------------------------- Save the merged count matrix ---------------------------##

write.table(count.matricies.merged, file= paste0(output.file.prefix,".merged.count.matrix.txt"), sep = "\t", quote = F, col.names = T, row.names = T)

##--------------------------- Regenerate basic intron metadata file and re-cluster ---------------------------##

intron.meta <- as.data.frame(str_split_fixed(count.matricies.merged$intron_junction, ":", 4), stringsAsFactors = FALSE )
names(intron.meta) <- c("chr","start","end","strand")

print(intron.meta[1:10,])

##Adjusting the meta-data taking strandedness into account
intron.meta = intron.meta %>% mutate(five.prime = ifelse(strand=="+", start, end), three.prime = ifelse(strand=="+", end, start))

#Five prime groups
fp.groups = intron.meta %>% group_indices(chr, five.prime)
intron.meta$clusterID_5p <- paste("clu_" ,fp.groups, sep="")

#Three prime groups
tp.groups = intron.meta %>% group_indices(chr, three.prime)
intron.meta$clusterID_3p <- paste("clu_" ,tp.groups, sep="")

print(intron.meta[1:10,])



write.table(intron.meta, file=paste0(output.file.prefix,".merged.intron.metadata.txt"), sep = "\t", quote = F, col.names = T, row.names = T)





