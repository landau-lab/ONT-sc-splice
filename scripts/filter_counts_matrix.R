## filter counts matrix and output filtered version of counts matrix and corresponding metadata 

args = commandArgs(TRUE)
path.to.counts = args[1]
path.to.metadata = args[2]
out.dir = args[3]
reads = args[4]

counts = read.table(path.to.counts)
#counts = read.table("/gpfs/commons/home/ahawkins/MDS_ONT_splice/leafcutter_longread/p1/full_matrix_no_3p_filter/strand_adjusted_matrix/strand_adjusted_counts.txt")

counts$junc_total = rowSums(counts)
counts.filt = counts[which(counts$junc_total >= reads),]

metadata = read.table(path.to.metadata)
#metadata = read.table("/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p2/2.Strand_Adjusted_Matrix/strand_adjusted_metadata.txt")
#counts = read.table("/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p2/3.Counts_Filtered_Matrix/counts.filtered.txt")

metadata$intron_junction_2 = paste(metadata$intron_junction, metadata$clusterID, sep=":")
rownames(metadata) = metadata$intron_junction_2
shared = intersect(metadata$intron_junction_2, rownames(counts.filt))
metadata.filt = metadata[shared,]

setwd(out.dir)
#setwd("/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/p2/3.Counts_Filtered_Matrix/")
write.table(counts.filt, "counts.filtered.txt")
write.table(metadata.filt, "metadata.filtered.txt", quote= FALSE, sep = "\t")
