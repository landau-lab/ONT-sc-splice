### merge outputs from parallelized permutation, combined patient
library(tidyverse)

args = commandArgs(TRUE)
path.to.outputs = args[1]
path.to.metadata = args[2]
output = args[3]


## test data
# path.to.outputs = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/MDS_P5_P6_2WT/5_read_threshold/diff_transcript_combined_merge_counts_ind_celltypes_2WT/NP/split_cluster_output"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/MDS_P5_P6_2WT/5_read_threshold/diff_transcript_combined_merge_counts_ind_celltypes_2WT/combined_metadata/combined_metadata.csv"

path.to.three.outputs = paste(path.to.outputs, "/alt_three_prime/", sep = "")
path.to.five.outputs = paste(path.to.outputs, "/alt_five_prime/", sep = "")

## read in all outputs 
setwd(path.to.three.outputs)
files = list.files(path.to.three.outputs)
three.output.list = lapply(files, function(x) read_csv(file=x))
three.output = do.call(rbind, three.output.list)

setwd(path.to.five.outputs)
files = list.files(path.to.five.outputs)
five.output.list = lapply(files, function(x) read.csv(file=x))
five.output = do.call(rbind, five.output.list)
message("output files loaded")
# 
# five.output = five.output %>% select(-five.wt.cluster.cov, -five.mut.cluster.cov)
# five.cluster.cov = five.output %>% group_by(three_prime_ID) %>% summarise(five.mut.cluster.cov = sum(obs.mut), five.wt.cluster.cov = sum(obs.wt))
# five.output = left_join(five.output, five.cluster.cov)

## load in metadata
metadata = read.csv(path.to.metadata)
metadata$intron_junction = paste(metadata$chr, metadata$start, metadata$end, sep = ":")
metadata = metadata %>% select(-patient)
metadata = distinct(metadata)
metadata$intron_junction = paste(metadata$intron_junction, metadata$strand, sep = ":")
metadata = metadata %>% select( - three_prime_ID, -five_prime_ID, -three_prime, -five_prime)

## merge output with metadata 
three.merge.output = inner_join(three.output, metadata)
five.merge.output = inner_join(five.output, metadata)
message("Output merged with metadata")

## calculate PSI for wt and mut 

three.merge.output$wt.psi = (three.merge.output$obs.wt/three.merge.output$three.wt.cluster.cov)*100
three.merge.output$mut.psi = (three.merge.output$obs.mut/three.merge.output$three.mut.cluster.cov)*100
three.merge.output$dPSI = three.merge.output$mut.psi - three.merge.output$wt.psi

five.merge.output$wt.psi = (five.merge.output$obs.wt/five.merge.output$five.wt.cluster.cov)*100
five.merge.output$mut.psi = (five.merge.output$obs.mut/five.merge.output$five.mut.cluster.cov)*100
five.merge.output$dPSI = five.merge.output$mut.psi - five.merge.output$wt.psi

message("PSI calculated")

## create heatmaps 
# sig.three = three.prime.data %>% filter(pvalue < 0.05)
# cryptic.five = five.prime.data %>% filter((startVerdict == "cryptic_fiveprime" | endVerdict == "cryptic_fiveprime") & pvalue < 0.05)
# 
# heat.three = sig.three[, grep("dPSI", colnames(sig.three))]
# heat.three$labels = paste(sig.three$gene, sig.three$intron_junction, sep = ":")
# heat.three = distinct(heat.three)
# labels = heat.three$labels
# 
# heat.three = as.matrix(heat.three[,1:3])
# rownames(heat.three) = labels
# 
# heat.three = replace(heat.three, is.na(heat.three), 0)
# 
# pheatmap(heat.three, fontsize_row = 5, main = "MDS Patients, All Junctions with p < 0.05 1000 perm")

## write output 
# setwd("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/4.logOR.combined.analysis/1000_within_cell_type/MDS_P1_P2_P3")
setwd(output)
write.table(three.merge.output, "logOR_within_cell_type_ALT_3P_Junctions.txt")
write.table(five.merge.output, "logOR_within_cell_type_ALT_5P_Junctions.txt")
message("DONE!")
