## merge combined cell type output for individual patients

library(tidyverse)

args = commandArgs(TRUE)
path.to.outputs = args[1]
path.to.metadata = args[2]
path.to.genotype = args[3]
output = args[4]


## test data
# path.to.outputs = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p2_biotin/ONT_pipeline/output_files/diff_transcript_output/split_cluster_celltype_output"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/p2_biotin/ONT_pipeline/output_files/strand_adjusted_metadata/strand_adjusted_metadata.csv"
# path.to.genotype = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p2.genotype.info.with.cell.type.txt"

path.to.three.outputs = paste(path.to.outputs, "/alt_three_prime/", sep = "")
path.to.five.outputs = paste(path.to.outputs, "/alt_five_prime/", sep = "")

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

## merge output with metadata 
three.merge.output = inner_join(three.output, metadata, by = "intron_junction")
five.merge.output = inner_join(five.output, metadata, by = "intron_junction")
message("Output merged with metadata")

## calculate psi for each cell type 

## get list of all cell types 
genotype = read.table(path.to.genotype)
genotype$Cell.Assignment = as.character(genotype$Cell.Assignment)
genotype[grep("EP_", genotype$Cell.Assignment), "Cell.Assignment"] = "EP"
cell.types = genotype$Cell.Assignment
cell.types = as.character(unique(cell.types))
cell.types[is.na(cell.types)] = "None"

three.cell.type.psi = data.frame(intron_junction = three.merge.output$intron_junction)
five.cell.type.psi = data.frame(intron_junction = five.merge.output$intron_junction)

for (cell in cell.types){
  wt.psi = (three.merge.output[, grep(paste("^", paste(cell, "obs.wt", sep = "_"), sep = ""), colnames(three.merge.output))]/three.merge.output[, grep(paste("^", paste(cell, "three.wt.cluster.cov", sep = "_"), sep = ""), colnames(three.merge.output))])*100
  mut.psi = (three.merge.output[, grep(paste("^", paste(cell, "obs.mut", sep = "_"), sep = ""), colnames(three.merge.output))]/three.merge.output[, grep(paste("^", paste(cell, "three.mut.cluster.cov", sep = "_"), sep = ""), colnames(three.merge.output))])*100
  if (length(wt.psi)[1] > 0 & length(mut.psi)[1] > 0){
    dPSI = mut.psi - wt.psi
    cell.data = data.frame(wt.psi = wt.psi, mut.psi = mut.psi, dPSI = dPSI)
    colnames(cell.data) = lapply(unlist(colnames(cell.data)), function(x) paste(cell, x, sep = "_"))
    cell.data$intron_junction = three.merge.output$intron_junction
    three.cell.type.psi = left_join(three.cell.type.psi, cell.data, by = "intron_junction")
  }
}

for (cell in cell.types){
  wt.psi = (five.merge.output[, grep(paste("^", paste(cell, "obs.wt", sep = "_"), sep = ""), colnames(five.merge.output))]/five.merge.output[, grep(paste("^",paste(cell, "five.wt.cluster.cov", sep = "_"), sep = ""), colnames(five.merge.output))])*100
  mut.psi = (five.merge.output[, grep(paste("^", paste(cell, "obs.mut", sep = "_"), sep = ""), colnames(five.merge.output))]/five.merge.output[, grep(paste("^", paste(cell, "five.mut.cluster.cov", sep = "_"), sep = ""), colnames(five.merge.output))])*100
  if (length(wt.psi)[1] > 0 & length(mut.psi)[1] > 0){
    dPSI = mut.psi - wt.psi
    cell.data = data.frame(wt.psi = wt.psi, mut.psi = mut.psi, dPSI = dPSI)
    colnames(cell.data) = lapply(unlist(colnames(cell.data)), function(x) paste(cell, x, sep = "_"))
    cell.data$intron_junction = five.merge.output$intron_junction
    five.cell.type.psi = left_join(five.cell.type.psi, cell.data, by = "intron_junction")
  }
}

message("PSI calculated")

## merge PSI information with final output

three.merge.output = inner_join(three.merge.output, three.cell.type.psi, by = "intron_junction")
five.merge.output = inner_join(five.merge.output, five.cell.type.psi, by = "intron_junction")

### add in annotation columns 
## add Final_Verdict column 
three.merge.output$Final_Verdict = NA
three.merge.output[which(three.merge.output$startClass == "main" & three.merge.output$endClass =="main"),"Final_Verdict"] = "Canonical"
three.merge.output[which(((three.merge.output$startClass == "main" & three.merge.output$endClass =="not_main_3_prime") | (three.merge.output$startClass == "not_main_3_prime" & three.merge.output$endClass =="main")) &
                   (three.merge.output$threep_distance > (-100) & three.merge.output$threep_distance < 0)),"Final_Verdict"] = "Cryptic_threeprime"
three.merge.output[which(((three.merge.output$startClass == "main" & three.merge.output$endClass =="not_main_5_prime") | (three.merge.output$startClass == "not_main_5_prime" & three.merge.output$endClass =="main")) &
                   (three.merge.output$fivep_distance > (-100) & three.merge.output$fivep_distance < 0)),"Final_Verdict"] = "Cryptic_fiveprime"
three.merge.output[which(three.merge.output$startClass == "not_main_3_prime" & three.merge.output$endClass =="not_main_3_prime"),"Final_Verdict"] = "cryptic_unanchored"
three.merge.output[which(three.merge.output$startClass == "not_main_5_prime" & three.merge.output$endClass =="not_main_5_prime"),"Final_Verdict"] = "cryptic_unanchored"
table(three.merge.output$Final_Verdict)

five.merge.output$merge.output_Verdict = NA
five.merge.output[which(five.merge.output$startClass == "main" & five.merge.output$endClass =="main"),"Final_Verdict"] = "Canonical"
five.merge.output[which(((five.merge.output$startClass == "main" & five.merge.output$endClass =="not_main_3_prime") | (five.merge.output$startClass == "not_main_3_prime" & five.merge.output$endClass =="main")) &
                  (five.merge.output$threep_distance > (-100) & five.merge.output$threep_distance < 0)),"Final_Verdict"] = "Cryptic_threeprime"
five.merge.output[which(((five.merge.output$startClass == "main" & five.merge.output$endClass =="not_main_5_prime") | (five.merge.output$startClass == "not_main_5_prime" & five.merge.output$endClass =="main")) &
                  (five.merge.output$fivep_distance > (-100) & five.merge.output$fivep_distance < 0)),"Final_Verdict"] = "Cryptic_fiveprime"
five.merge.output[which(five.merge.output$startClass == "not_main_3_prime" & five.merge.output$endClass =="not_main_3_prime"),"Final_Verdict"] = "cryptic_unanchored"
five.merge.output[which(five.merge.output$startClass == "not_main_5_prime" & five.merge.output$endClass =="not_main_5_prime"),"Final_Verdict"] = "cryptic_unanchored"
table(five.merge.output$Final_Verdict)

## add in threshold columns 
# three.merge.output$dPSI_threshold_0 = "no"
# five.merge.output$dPSI_threshold_0 = "no"
# three.merge.output$dPSI_threshold_5 = "no"
# five.merge.output$dPSI_threshold_5 = "no"
# three.merge.output$pvalue_threshold = "no"
# five.merge.output$pvalue_threshold = "no"
# three.merge.output[which(three.merge.output$HSPC_dPSI > 0 | three.merge.output$IMP_dPSI > 0 | three.merge.output$MEP_dPSI > 0 |
#                            three.merge.output$HSPC_EP_dPSI > 0 | three.merge.output$MkP_dPSI > 0 | three.merge.output$EP_1_dPSI > 0 |
#                            three.merge.output$NP_dPSI > 0 | three.merge.output$CC_dPSI > 0 | three.merge.output$IMP_MEP_dPSI > 0),]$dPSI_threshold_0 = "yes"
# five.merge.output[which(five.merge.output$HSPC_dPSI > 0 | five.merge.output$IMP_dPSI > 0 | five.merge.output$MEP_dPSI > 0 |
#                            five.merge.output$HSPC_EP_dPSI > 0 | five.merge.output$MkP_dPSI > 0 | five.merge.output$EP_1_dPSI > 0 |
#                            five.merge.output$NP_dPSI > 0 | five.merge.output$CC_dPSI > 0 | five.merge.output$IMP_MEP_dPSI > 0),]$dPSI_threshold_0 = "yes"
# three.merge.output[which(three.merge.output$HSPC_dPSI >= 5 | three.merge.output$IMP_dPSI >= 5 | three.merge.output$MEP_dPSI >= 5 |
#                            three.merge.output$HSPC_EP_dPSI >= 5 | three.merge.output$MkP_dPSI >= 5 | three.merge.output$EP_1_dPSI >= 5 |
#                            three.merge.output$NP_dPSI >= 5 | three.merge.output$CC_dPSI >= 5 | three.merge.output$IMP_MEP_dPSI >= 5),]$dPSI_threshold_5 = "yes"
# five.merge.output[which(five.merge.output$HSPC_dPSI >= 5 | five.merge.output$IMP_dPSI >= 5 | five.merge.output$MEP_dPSI >= 5 |
#                           five.merge.output$HSPC_EP_dPSI >= 5 | five.merge.output$MkP_dPSI >= 5 | five.merge.output$EP_1_dPSI >= 5 |
#                           five.merge.output$NP_dPSI >= 5 | five.merge.output$CC_dPSI >= 5 | five.merge.output$IMP_MEP_dPSI >= 5),]$dPSI_threshold_5 = "yes"
# three.merge.output[which(three.merge.output$HSPC_pvalue < 0.05 | three.merge.output$IMP_pvalue < 0.05 | three.merge.output$MEP_pvalue < 0.05 |
#                            three.merge.output$HSPC_EP_pvalue < 0.05 | three.merge.output$MkP_pvalue < 0.05 | three.merge.output$EP_1_pvalue < 0.05 |
#                            three.merge.output$NP_pvalue < 0.05 | three.merge.output$CC_pvalue < 0.05 | three.merge.output$IMP_MEP_pvalue < 0.05),]$pvalue_threshold = "yes"
# five.merge.output[which(five.merge.output$HSPC_pvalue < 0.05 | five.merge.output$IMP_pvalue < 0.05 | five.merge.output$MEP_pvalue < 0.05 |
#                           five.merge.output$HSPC_EP_pvalue < 0.05 | five.merge.output$MkP_pvalue < 0.05 | five.merge.output$EP_1_pvalue < 0.05 |
#                           five.merge.output$NP_pvalue < 0.05 | five.merge.output$CC_pvalue < 0.05 | five.merge.output$IMP_MEP_pvalue < 0.05),]$pvalue_threshold = "yes"
# 
# 
# 
# 
# ## add in dPSI and pvalue double threshold columns
# three.merge.output$dPSI_0_pvalue_threshold = "no"
# five.merge.output$dPSI_0_pvalue_threshold = "no"
# three.merge.output$dPSI_5_pvalue_threshold = "no"
# five.merge.output$dPSI_5_pvalue_threshold = "no"
# three.merge.output[which((three.merge.output$HSPC_dPSI > 0 & three.merge.output$HSPC_pvalue < 0.05) | (three.merge.output$IMP_dPSI > 0 & three.merge.output$IMP_pvalue < 0.05) | (three.merge.output$MEP_dPSI > 0 & three.merge.output$MEP_pvalue < 0.05)|
#                            (three.merge.output$HSPC_EP_dPSI > 0 & three.merge.output$HSPC_EP_pvalue < 0.05) | (three.merge.output$MkP_dPSI > 0 & three.merge.output$MkP_pvalue < 0.05) | (three.merge.output$EP_1_dPSI > 0 & three.merge.output$EP_1_pvalue < 0.05) |
#                            (three.merge.output$NP_dPSI > 0 & three.merge.output$NP_pvalue < 0.05) | (three.merge.output$CC_dPSI > 0 & three.merge.output$CC_pvalue < 0.05) | (three.merge.output$IMP_MEP_dPSI > 0 & three.merge.output$IMP_MEP_pvalue < 0.05)),]$dPSI_0_pvalue_threshold = "yes"
# five.merge.output[which((five.merge.output$HSPC_dPSI > 0 & five.merge.output$HSPC_pvalue < 0.05) | (five.merge.output$IMP_dPSI > 0 & five.merge.output$IMP_pvalue < 0.05) | (five.merge.output$MEP_dPSI > 0 & five.merge.output$MEP_pvalue < 0.05)|
#                            (five.merge.output$HSPC_EP_dPSI > 0 & five.merge.output$HSPC_EP_pvalue < 0.05) | (five.merge.output$MkP_dPSI > 0 & five.merge.output$MkP_pvalue < 0.05) | (five.merge.output$EP_1_dPSI > 0 & five.merge.output$EP_1_pvalue < 0.05) |
#                            (five.merge.output$NP_dPSI > 0 & five.merge.output$NP_pvalue < 0.05) | (five.merge.output$CC_dPSI > 0 & five.merge.output$CC_pvalue < 0.05) | (five.merge.output$IMP_MEP_dPSI > 0 & five.merge.output$IMP_MEP_pvalue < 0.05)),]$dPSI_0_pvalue_threshold = "yes"
# three.merge.output[which((three.merge.output$HSPC_dPSI >= 5 & three.merge.output$HSPC_pvalue < 0.05) | (three.merge.output$IMP_dPSI >= 5 & three.merge.output$IMP_pvalue < 0.05) | (three.merge.output$MEP_dPSI >= 5 & three.merge.output$MEP_pvalue < 0.05)|
#                            (three.merge.output$HSPC_EP_dPSI >= 5 & three.merge.output$HSPC_EP_pvalue < 0.05) | (three.merge.output$MkP_dPSI >= 5 & three.merge.output$MkP_pvalue < 0.05) | (three.merge.output$EP_1_dPSI >= 5 & three.merge.output$EP_1_pvalue < 0.05) |
#                            (three.merge.output$NP_dPSI >= 5 & three.merge.output$NP_pvalue < 0.05) | (three.merge.output$CC_dPSI >= 5 & three.merge.output$CC_pvalue < 0.05) | (three.merge.output$IMP_MEP_dPSI >= 5 & three.merge.output$IMP_MEP_pvalue < 0.05)),]$dPSI_5_pvalue_threshold = "yes"
# five.merge.output[which((five.merge.output$HSPC_dPSI >= 5 & five.merge.output$HSPC_pvalue < 0.05) | (five.merge.output$IMP_dPSI >= 5 & five.merge.output$IMP_pvalue < 0.05) | (five.merge.output$MEP_dPSI >= 5 & five.merge.output$MEP_pvalue < 0.05)|
#                           (five.merge.output$HSPC_EP_dPSI >= 5 & five.merge.output$HSPC_EP_pvalue < 0.05) | (five.merge.output$MkP_dPSI >= 5 & five.merge.output$MkP_pvalue < 0.05) | (five.merge.output$EP_1_dPSI >= 5 & five.merge.output$EP_1_pvalue < 0.05) |
#                           (five.merge.output$NP_dPSI >= 5 & five.merge.output$NP_pvalue < 0.05) | (five.merge.output$CC_dPSI >= 5 & five.merge.output$CC_pvalue < 0.05) | (five.merge.output$IMP_MEP_dPSI >= 5 & five.merge.output$IMP_MEP_pvalue < 0.05)),]$dPSI_5_pvalue_threshold = "yes"

## filter for cryptic junctions only 
three.final.cryptic = three.merge.output %>% filter(Final_Verdict == "Cryptic_threeprime")
five.final.cryptic = five.merge.output %>% filter(Final_Verdict == "Cryptic_fiveprime")

## write output 
setwd(output)
write.table(three.merge.output, "all_celltypes_DTU_Alt_3P_withPSI.txt")
write.table(five.merge.output, "all_celltypes_DTU_Alt_5P_withPSI.txt")
write.table(three.final.cryptic, "all_celltypes_DTU_CRYPTIC_3P_withPSI.txt")
write.table(five.final.cryptic, "all_celltypes_DTU_CRYPTIC_5P_withPSI.txt")
