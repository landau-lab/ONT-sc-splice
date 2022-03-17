###############################
#### this script generates the input matrices to do differential transcript usage ######
# input is full counts matrix and metadata
# output is strand adjusted counts matrix and metadata
# argument 1 is the full counts matrix from leafcutter
# argument 2 is the full annotated matrix from leafcutter (with additional annotations from Lloyd's annotation script)
# argument 3 is the directory to output the files 
###############################

library(tidyverse)

#load directories to input files 
args = commandArgs(TRUE)
path.to.metadata = args[1]
output.dir = args[2]

##test data
# path.to.metadata = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/leafcutter_outputs/MDS_P1_output/MDS_P1_all.introns.info.w.primary.annotations.txt"
# output.dir = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p1_biotin/strand_adjusted_outputs"

#load in metadata 
all_introns_meta = read.csv(path.to.metadata, row.names = 1)
print("Adjusting for strandedness...")

#adjust for strandedness
all_introns_meta$fivep_class = 'NA'
all_introns_meta$threep_class = 'NA'
all_introns_meta = all_introns_meta %>% mutate(fivep_class = ifelse(strand_1=="+", startClass, endClass),
                                               threep_class = ifelse(strand_1 == "+", endClass, startClass))
print("Taking strandedness into account")

#Taking strandedness into account 
all_introns_meta$fivep_distance <- NA 
all_introns_meta$threep_distance <- NA
all_introns_meta[which(all_introns_meta$strand_1 == "+"), ]$fivep_distance <- all_introns_meta[which(all_introns_meta$strand_1 == "+"), ]$startDistance
all_introns_meta[which(all_introns_meta$strand_1 == "+"), ]$threep_distance <- all_introns_meta[which(all_introns_meta$strand_1 == "+"), ]$endDistance
all_introns_meta[which(all_introns_meta$strand_1 == "-"),]$fivep_distance <- (all_introns_meta[which(all_introns_meta$strand_1 == "-"), ]$endDistance)*-1
all_introns_meta[which(all_introns_meta$strand_1 == "-"),]$threep_distance <- (all_introns_meta[which(all_introns_meta$strand_1 == "-"), ]$startDistance)*-1
print("Strand adjustment done")
print(dim(all_introns_meta))
#subset metadata by those that have strand information and filter counts information for only those that have strand information 

all_introns_meta$intron_junction = paste(all_introns_meta$chr, all_introns_meta$start, all_introns_meta$end, all_introns_meta$strand_1, sep = ":")

all_introns_meta = all_introns_meta %>% filter(strand_1 != 'NA')

print("Writing output files...")
setwd(output.dir)
write_csv(all_introns_meta, "strand_adjusted_metadata.csv")
print("Done")
