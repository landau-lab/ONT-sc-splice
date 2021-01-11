### merge outputs from parallelized permutation, combined patient
library(tidyverse)

args = commandArgs(TRUE)
path.to.outputs = args[1]
path.to.metadata = args[2]
patient.names = args[3]
output = args[4]


## test data
# path.to.outputs = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output_remove_celltypes/NP/split_cluster_output"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/9.Splice_pipeline_example_run/CH_combined/diff_transcript_combined_output/combined_metadata/combined_metadata.csv"

path.to.three.outputs = paste(path.to.outputs, "/alt_three_prime/", sep = "")
path.to.five.outputs = paste(path.to.outputs, "/alt_five_prime/", sep = "")
patient.names = unlist(strsplit(patient.names, split = ","))

## read in all outputs 
setwd(path.to.three.outputs)
files = list.files(path.to.three.outputs)
three.output.list = lapply(files, function(x) read.csv(file=x))
three.output = do.call(rbind, three.output.list)

setwd(path.to.five.outputs)
files = list.files(path.to.five.outputs)
five.output.list = lapply(files, function(x) read.csv(file=x))
five.output = do.call(rbind, five.output.list)
message("output files loaded")

## load in metadata
metadata = read.csv(path.to.metadata)
metadata$intron_junction = paste(metadata$chr, metadata$start, metadata$end, sep = ":")
metadata = metadata %>% select(-patient)
metadata = distinct(metadata)
metadata$intron_junction = paste(metadata$intron_junction, metadata$strand, sep = ":")

## merge output with metadata 
three.merge.output = inner_join(three.output, metadata, by = "intron_junction")
five.merge.output = inner_join(five.output, metadata, by = "intron_junction")
message("Output merged with metadata")

## calculate PSI for wt and mut for each patient 
three.split.out = map(set_names(patient.names),~select(three.merge.output,starts_with(.x)))
for (patient in patient.names){
  colnames(three.split.out[[patient]]) = gsub(paste(patient, "_", sep = ""), "", colnames(three.split.out[[patient]]))
  cluster.cov.wt = three.split.out[[patient]] %>% group_by(five_prime_ID) %>% 
    summarise("wt.cluster.cov" = sum(obs.wt))
  cluster.cov.mut = three.split.out[[patient]] %>% group_by(five_prime_ID) %>% 
    summarise("mut.cluster.cov" = sum(obs.mut))
  three.split.out[[patient]] = left_join(three.split.out[[patient]], cluster.cov.wt, by = "five_prime_ID")
  three.split.out[[patient]] = left_join(three.split.out[[patient]], cluster.cov.mut, by = "five_prime_ID")
  three.split.out[[patient]]$wt.psi = (three.split.out[[patient]]$obs.wt/three.split.out[[patient]]$wt.cluster.cov)*100
  three.split.out[[patient]]$mut.psi = (three.split.out[[patient]]$obs.mut/three.split.out[[patient]]$mut.cluster.cov)*100
  three.split.out[[patient]]$dPSI = three.split.out[[patient]]$mut.psi - three.split.out[[patient]]$wt.psi
  colnames(three.split.out[[patient]])= lapply(colnames(three.split.out[[patient]]), function(x) paste(patient,x,sep="_"))
}

five.split.out = map(set_names(patient.names),~select(five.merge.output,starts_with(.x)))
for (patient in patient.names){
  colnames(five.split.out[[patient]]) = gsub(paste(patient, "_", sep = ""), "", colnames(five.split.out[[patient]]))
  cluster.cov.wt = five.split.out[[patient]] %>% group_by(three_prime_ID) %>% 
    summarise("wt.cluster.cov" = sum(obs.wt))
  cluster.cov.mut = five.split.out[[patient]] %>% group_by(three_prime_ID) %>% 
    summarise("mut.cluster.cov" = sum(obs.mut))
  five.split.out[[patient]] = left_join(five.split.out[[patient]], cluster.cov.wt)
  five.split.out[[patient]] = left_join(five.split.out[[patient]], cluster.cov.mut)
  five.split.out[[patient]]$wt.psi = (five.split.out[[patient]]$obs.wt/five.split.out[[patient]]$wt.cluster.cov)*100
  five.split.out[[patient]]$mut.psi = (five.split.out[[patient]]$obs.mut/five.split.out[[patient]]$mut.cluster.cov)*100
  five.split.out[[patient]]$dPSI = five.split.out[[patient]]$mut.psi - five.split.out[[patient]]$wt.psi
  colnames(five.split.out[[patient]])= lapply(colnames(five.split.out[[patient]]), function(x) paste(patient,x,sep="_"))
}

## merge back together into one data frame 
three.out = bind_cols(three.split.out)
three.prime.data = inner_join(three.merge.output, three.out)


five.out = bind_cols(five.split.out)
five.prime.data = inner_join(five.merge.output, five.out)

## calculate wieghted PSI 
## now calculated weighted average of PSI 
# library(matrixStats)
# 
# 
# three.total.junction.coverage = three.prime.data[,grep("total.reads.per.junction", colnames(three.prime.data))]
# five.total.junction.coverage = five.prime.data[,grep("total.reads.per.junction", colnames(five.prime.data))]
# three.weights = three.total.junction.coverage/rowSums(three.total.junction.coverage)
# five.weights = five.total.junction.coverage/rowSums(five.total.junction.coverage)
# 
# 
# 
# three.wt.psi = three.prime.data[,grep("wt.psi", colnames(three.prime.data))]
# three.mut.psi = three.prime.data[,grep("mut.psi", colnames(three.prime.data))]
# i = 1
# three.prime.data$weight.wt.PSI = NA
# three.prime.data$weight.mut.PSI = NA
# while (i <= dim(three.wt.psi)[1]){
#   three.prime.data[i, "weight.wt.PSI"] = rowWeightedMeans(as.matrix(three.wt.psi[i,]), as.numeric(three.weights[i,]))
#   three.prime.data[i, "weight.mut.PSI"] = rowWeightedMeans(as.matrix(three.mut.psi[i,]), as.numeric(three.weights[i,]))
#   i = i +1
# }
# 
# three.prime.data$weight.dPSI = three.prime.data$weight.mut.PSI - three.prime.data$weight.wt.PSI
# 
# five.wt.psi = five.prime.data[,grep("wt.psi", colnames(five.prime.data))]
# five.mut.psi = five.prime.data[,grep("mut.psi", colnames(five.prime.data))]
# 
# i = 1
# five.prime.data$weight.wt.PSI = NA
# five.prime.data$weight.mut.PSI = NA
# while (i <= dim(five.wt.psi)[1]){
#   five.prime.data[i, "weight.wt.PSI"] = rowWeightedMeans(as.matrix(five.wt.psi[i,]), as.numeric(five.weights[i,]))
#   five.prime.data[i, "weight.mut.PSI"] = rowWeightedMeans(as.matrix(five.mut.psi[i,]), as.numeric(five.weights[i,]))
#   i = i +1
# }
# 
# five.prime.data$weight.dPSI = five.prime.data$weight.mut.PSI - five.prime.data$weight.wt.PSI
# message("PSI calculated")

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
write.table(three.prime.data, "combined_patient_logOR_within_cell_type_ALT_3P_Junctions.txt")
write.table(five.prime.data, "combined_patient_logOR_within_cell_type_ALT_5P_Junctions.txt")
message("DONE!")
