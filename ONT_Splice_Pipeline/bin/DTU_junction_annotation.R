## add in annotations for categories for each junction 
## columns being added in are as follows:
## Final_Verdict - Canonical, Cryptic Threeprime, Cryptic Fiveprime, or Cryptic Unanchored 
## dPSI_threshold_0 - Cryptic threeprime/ fiveprime + dPSI > 0 
## dPSI_threshold_5 - Cryptic threeprime/ fiveprime + dPSI >= 5
## pvalue_threshold - Cryptic threeprime/ fiveprime + pvalue < 0.05
## dPSI_0_pvalue_threshold - Cryptic threeprime/ fiveprime + dPSI > 0, pvalue < 0.05
## dPSI_5_pvalue_threshold - Cryptic threeprime/ fiveprime + dPSI >= 5, pvalue < 0.05 

library(tidyverse)

args = commandArgs(TRUE)
path.to.three.data = args[1]
path.to.five.data = args[2]
output = args[3]

## test data
# path.to.three.data = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/3.combined_patients/CH_combined/all_celltypes/combined_patient_logOR_within_cell_type_ALT_3P_Junctions.txt"
# path.to.five.data = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/3.combined_patients/CH_combined/all_celltypes/combined_patient_logOR_within_cell_type_ALT_5P_Junctions.txt"

three.data = read.table(path.to.three.data)
five.data = read.table(path.to.five.data)

## add Final_Verdict column 
three.data$Final_Verdict = NA
three.data[which(three.data$startClass == "main" & three.data$endClass =="main"),"Final_Verdict"] = "Canonical"
three.data[which(((three.data$startClass == "main" & three.data$endClass =="not_main_3_prime") | (three.data$startClass == "not_main_3_prime" & three.data$endClass =="main")) &
             (three.data$threep_distance > (-100) & three.data$threep_distance < 0)),"Final_Verdict"] = "Cryptic_threeprime"
three.data[which(((three.data$startClass == "main" & three.data$endClass =="not_main_5_prime") | (three.data$startClass == "not_main_5_prime" & three.data$endClass =="main")) &
                   (three.data$fivep_distance > (-100) & three.data$fivep_distance < 0)),"Final_Verdict"] = "Cryptic_fiveprime"
three.data[which(three.data$startClass == "not_main_3_prime" & three.data$endClass =="not_main_3_prime"),"Final_Verdict"] = "cryptic_unanchored"
three.data[which(three.data$startClass == "not_main_5_prime" & three.data$endClass =="not_main_5_prime"),"Final_Verdict"] = "cryptic_unanchored"
table(three.data$Final_Verdict)

five.data$data_Verdict = NA
five.data[which(five.data$startClass == "main" & five.data$endClass =="main"),"Final_Verdict"] = "Canonical"
five.data[which(((five.data$startClass == "main" & five.data$endClass =="not_main_3_prime") | (five.data$startClass == "not_main_3_prime" & five.data$endClass =="main")) &
                   (five.data$threep_distance > (-100) & five.data$threep_distance < 0)),"Final_Verdict"] = "Cryptic_threeprime"
five.data[which(((five.data$startClass == "main" & five.data$endClass =="not_main_5_prime") | (five.data$startClass == "not_main_5_prime" & five.data$endClass =="main")) &
                   (five.data$fivep_distance > (-100) & five.data$fivep_distance < 0)),"Final_Verdict"] = "Cryptic_fiveprime"
five.data[which(five.data$startClass == "not_main_3_prime" & five.data$endClass =="not_main_3_prime"),"Final_Verdict"] = "cryptic_unanchored"
five.data[which(five.data$startClass == "not_main_5_prime" & five.data$endClass =="not_main_5_prime"),"Final_Verdict"] = "cryptic_unanchored"
table(five.data$Final_Verdict)

## add in dPSI and pvalue single threshold columns 
three.data$dPSI_threshold_0 = NA
five.data$dPSI_threshold_0 = NA
three.data$dPSI_threshold_5 = NA
five.data$dPSI_threshold_5 = NA
three.data$pvalue_threshold = NA
five.data$pvalue_threshold = NA
three.data = three.data %>% mutate(dPSI_threshold_0 = ifelse(dPSI > 0, "yes", "no"),
                            dPSI_threshold_5 = ifelse(dPSI >= 5, "yes", "no"),
                            pvalue_threshold = ifelse(pvalue < 0.05, "yes", "no"))
five.data = five.data %>% mutate(dPSI_threshold_0 = ifelse(dPSI > 0, "yes", "no"),
                                   dPSI_threshold_5 = ifelse(dPSI >= 5, "yes", "no"),
                                 pvalue_threshold = ifelse(pvalue < 0.05, "yes", "no"))


## add in dPSI and pvalue double threshold columns
three.data$dPSI_0_pvalue_threshold = "no"
five.data$dPSI_0_pvalue_threshold = "no"
three.data$dPSI_5_pvalue_threshold = "no"
five.data$dPSI_5_pvalue_threshold = "no"
three.data[which(three.data$pvalue < 0.05 & three.data$dPSI > 0),]$dPSI_0_pvalue_threshold = "yes"
three.data[which(three.data$pvalue < 0.05 & three.data$dPSI >= 5),]$dPSI_5_pvalue_threshold = "yes"
five.data[which(five.data$pvalue < 0.05 & five.data$dPSI > 0),]$dPSI_0_pvalue_threshold = "yes"
five.data[which(five.data$pvalue < 0.05 & five.data$dPSI >= 5),]$dPSI_5_pvalue_threshold = "yes"

## filter for only cryptic sites 
three.data.cryptic = three.data %>% filter(Final_Verdict == "Cryptic_threeprime")
five.data.cryptic = five.data %>% filter(Final_Verdict == "Cryptic_fiveprime")

## make a plot to look at cryptic sites falling into different thresholds 
dpsi_0_threshold = as.data.frame(table(three.data.cryptic$dPSI_threshold_0))
dpsi_5_threshold = as.data.frame(table(three.data.cryptic$dPSI_threshold_5))
pvalue_threshold = as.data.frame(table(three.data.cryptic$pvalue_threshold))
dpsi_0_pvalue_threshold = as.data.frame(table(three.data.cryptic$dPSI_0_pvalue_threshold))
dpsi_5_pvalue_threshold = as.data.frame(table(three.data.cryptic$dPSI_5_pvalue_threshold))

three.df.list = list(dpsi_0_threshold, dpsi_5_threshold, pvalue_threshold, dpsi_0_pvalue_threshold, dpsi_5_pvalue_threshold)
names(three.df.list) = c("dPSI_0", "dPSI_5", "pvalue_0.05", "dPSI_0_pvalue_0.05", "dPSI_5_pvalue_0.05" )

three.df = bind_rows(three.df.list, .id = "Category")


dpsi_0_threshold = as.data.frame(table(five.data.cryptic$dPSI_threshold_0))
dpsi_5_threshold = as.data.frame(table(five.data.cryptic$dPSI_threshold_5))
pvalue_threshold = as.data.frame(table(five.data.cryptic$pvalue_threshold))
dpsi_0_pvalue_threshold = as.data.frame(table(five.data.cryptic$dPSI_0_pvalue_threshold))
dpsi_5_pvalue_threshold = as.data.frame(table(five.data.cryptic$dPSI_5_pvalue_threshold))

five.df.list = list(dpsi_0_threshold, dpsi_5_threshold, pvalue_threshold, dpsi_0_pvalue_threshold, dpsi_5_pvalue_threshold)
names(five.df.list) = c("dPSI_0", "dPSI_5", "pvalue_0.05", "dPSI_0_pvalue_0.05", "dPSI_5_pvalue_0.05" )

five.df = bind_rows(five.df.list, .id = "Category")


## save plots and files 
library(forcats)

setwd(output)
three.df %>% mutate(Var1 = fct_reorder(Var1, Category)) %>% 
  ggplot(aes(x = Category, y = Freq, fill = Var1)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  ggtitle("MDS_P1 Threeprime Cryptic Junctions") + theme(plot.title = element_text(hjust = 0.5))

ggsave("Cryptic_threeprime_thresholds_distribution.pdf", width = 10, height = 10)

five.df %>% mutate(Var1 = fct_reorder(Var1, Category)) %>% 
  ggplot(aes(x = Category, y = Freq, fill = Var1)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  ggtitle("MDS_P1 Fiveprime Cryptic Junctions") + theme(plot.title = element_text(hjust = 0.5))

ggsave("Cryptic_fiveprime_thresholds_distribution.pdf", width = 10, height = 10)

write.table(three.data, "logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt")
write.table(five.data, "logOR_within_cell_type_ALT_5P_Junctions_with_threshold_info.txt")
write.table(three.data.cryptic, "logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info.txt")
write.table(five.data.cryptic, "logOR_within_cell_type_ONLY_CRYPTIC_5P_Junctions_with_threshold_info.txt")
