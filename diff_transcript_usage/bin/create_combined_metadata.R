### create combined metadata across all three patients 
library(tidyverse)

args = commandArgs(TRUE)
path.to.metadata = args[1]
output.dir = args[2]
patient.names = args[3]

patient.names = unlist(strsplit(patient.names, split = ","))

print(patient.names)

## test data
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/3.Annotated_metadata"
# patient.names = c("CH259", "CH305")


## load in metadata
setwd(path.to.metadata)
files = list.files(path.to.metadata)
metadata.list = lapply(files, function(x) read.csv(file = x))
names(metadata.list) = patient.names

## identify strand adjusted start and end site
for (patient in patient.names) {
  metadata.list[[patient]]$five_prime = "NA"
  metadata.list[[patient]]$three_prime = "NA"
  metadata.list[[patient]] = metadata.list[[patient]] %>% mutate(five_prime = ifelse(strand == "+", start, end),
                                      three_prime = ifelse(strand =="+", end,start))
}

## select only parts of metadata that we care about
for (patient in patient.names) {
  metadata.list[[patient]] = metadata.list[[patient]] %>% select(chr, start, end, strand, five_prime, three_prime, gene, endClass, startClass, verdict, fivep_distance, threep_distance)
}

#combine into one list of all intron junctions and create new three prime and five prime cluster IDs 
metadata = do.call(rbind, metadata.list)
metadata$five_prime_ID = metadata %>% group_by(chr, five_prime, strand) %>% group_indices()
metadata$three_prime_ID = metadata %>% group_by(chr, three_prime, strand) %>% group_indices()
metadata$five_prime_ID = gsub("^", "clu_", metadata$five_prime_ID)
metadata$three_prime_ID = gsub("^", "clu_", metadata$three_prime_ID)
metadata$patient = rownames(metadata)
metadata$patient = gsub("\\..*", "", metadata$patient)

## write output
setwd(output.dir)
write.csv(metadata, "combined_metadata.csv", quote = FALSE, row.names = FALSE)

#write.csv(metadata, "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/6.combined_metadata/MDS_combined_metadata.txt", quote = FALSE)
#metadata = read.csv("/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/6.combined_metadata/MDS_combined_metadata.txt", stringsAsFactors = FALSE, row.names = 1)

