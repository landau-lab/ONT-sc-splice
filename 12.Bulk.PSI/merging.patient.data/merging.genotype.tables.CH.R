
##------------------------------- Merging the CH259 + CH305 genotyping tables ------------------------------- ##

GT.1 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_259.genotype.info.with.cell.type.txt")

GT.2 <- read.table("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/3.Genotyping_info/ch_305.genotype.info.with.cell.type.txt")

GT.merged <- rbind(GT.1, GT.2)

print("merging CH genotype tables")
nrow(GT.1)
nrow(GT.2)
nrow(GT.merged)

write.table(GT.merged, file="", )
