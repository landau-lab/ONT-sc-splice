
##------------------------------- Merging the MDSP5 + MDSP6 genotyping tables ------------------------------- ##

GT.1 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p5_1.genotype.info.with.cell.type.txt")

GT.2 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p5_2.genotype.info.with.cell.type.txt")

GT.3 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/mds_p6.genotype.info.with.cell.type.txt")

GT.merged.temp <- rbind(GT.1, GT.2)

GT.merged <- rbind(GT.merged.temp, GT.3)

print("merging MDS genotype tables")
nrow(GT.1)
nrow(GT.2)
nrow(GT.3)
nrow(GT.merged)

write.table(GT.merged, file="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/merging.patient.data/MDSP5.MDSP6.merged.genotype.info.with.cell.type.txt", quote = F, row.names = T, col.names = T)
