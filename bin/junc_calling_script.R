require(tidyverse)
require(magrittr)
require(irlba)
library(optparse)
library(dplyr)
library(data.table)
library(stringr)


option_parser=OptionParser(
  usage="%prog [options] <name>_counts_sc_txt.gz path/to/output_folder output_prefix path/to/annotation_code/prefix_"
)

parsed_args <- parse_args(option_parser, positional_arguments = 4)

input_matrix <- parsed_args$args[1]
output_folder <- parsed_args$args[2]
output_prefix <- parsed_args$args[3]
annotation_code <- parsed_args$args[4]


cat("Results to be saved in:",output_folder, "\n")

output <- paste0(output_folder,"/",output_prefix)


dat = read_tsv(input_matrix)
print("Saving matrix as an R object")
saveRDS(dat, version=2, paste0(output,"_counts_sc.rda"))
print("Done")

#dat = readRDS("/gpfs/commons/home/pchamely/leafcutter_scripts/p1_output/p1_ont_v2_align_test/p1_align_smartseq_counts_sc.rda")

#junc_meta = dat %>% select("chrom", "strand", "start", "end")
junc_meta = dat[,1:4]
counts = dat[,5:ncol(dat)]

#Adding on the junction information to the counts matrix!
junc_meta$intron_junction <- paste(junc_meta$chrom, ":", junc_meta$start, ":", junc_meta$end,":",junc_meta$strand, sep = "")

#Save count matrix
perind_numbers.counts <- bind_cols(junc_meta[,c("intron_junction")], counts)
print("Writing out full count matrix table")
write.table(perind_numbers.counts, file= paste0(output ,"_perind_numbers.counts.txt") , row.names = F, col.names = T, quote = F, sep = " " )
print("Done")

#perind_numbers.counts[duplicated(perind_numers.counts$intron_junction),] %>% arrange(intron_junction)


## Junction Annotation:
print("Beginning Junction Annotation Process")

# annotation
#exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0(annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0(annotation_code,"_fiveprime.bed.gz")

cat("Using annotation at:", annotation_code,"\n")

##This is where it reads in the count matrix with all the junctions and cells
cat("Loading counts from", paste0(output ,"_perind_numers.counts.txt"), "\n")
#counts <- read.table(counts_file, check.names=FALSE)
counts <- perind_numbers.counts
counts <- as.data.frame(counts)
rownames(counts) <- counts$intron_junction
print("finished loading count matirx")

#For the test run
add_chr=function(chrs)
  if (!grepl("chr",chrs[1])) paste0("chr",chrs) else chrs

all.introns <- as.data.frame(str_split_fixed(rownames(counts), ":", 4), stringsAsFactors = FALSE )
names(all.introns) <- c("chr","start","end","strand_1")


##--------added in here------##
##Adjusting the meta-data taking strandedness into account
all.introns = all.introns %>% mutate(five_prime = ifelse(strand_1 =="+", start, end),
                                     three_prime = ifelse(strand_1 =="+", end, start))

#Five prime groups
fp_groups = all.introns %>% group_indices(chr, five_prime)
all.introns$clusterID_5p <- paste("clu_" ,fp_groups, sep="")

#Three prime groups
tp_groups = all.introns %>% group_indices(chr, three_prime)
all.introns$clusterID_3p <- paste("clu_" ,tp_groups, sep="")

##--------added in here------##


#all.introns$chr <- add_chr(all.introns$chr)
if(nrow(all.introns) == 0 ){
  stop("No intron regions found! Please check your input files.")
}

all.introns$start <- as.numeric(all.introns$start)
all.introns$end <- as.numeric(all.introns$end)
all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID_5p)
colnames(all.junctions)[4] <- "clusterID"

##Making an intron database
#paste("zcat < ",all_introns)
intron_db <- fread(all_introns, data.table = FALSE)
colnames(intron_db)[1:4]=c("chr","start","end","gene")
intron_db$end <- (intron_db$end-1)
intron_db$chr <- add_chr(intron_db$chr)
print("finished making intron database")

##Making the 3' database
threeprime_db <- fread(threeprime_file, data.table = FALSE)
colnames(threeprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
threeprime_db$end <- (threeprime_db$end-1)
threeprime_db$start <- (threeprime_db$start-1)
threeprime_db$chr <- add_chr(threeprime_db$chr)
print("finished making 3P database")

##Making the 5' database
#fiveprime_db <- fread(paste("zcat < ",fiveprime_file), data.table = FALSE)
fiveprime_db <- fread(fiveprime_file, data.table = FALSE)
colnames(fiveprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
fiveprime_db$chr <- add_chr(fiveprime_db$chr)
print("finished making 5' database")

##Adjusting the all_introns matrix
all.introns$start_match <- NA
all.introns$end_match <- NA

##try putting this before and if it works then see if you can create a new row with the adjusted start and end positions of the junctions!!!
for (chrom in unique(all.introns$chr)){

  all.introns[which(all.introns$chr == chrom),]$start_match <- sapply(all.introns[which(all.introns$chr == chrom),]$start, function(x) fiveprime_db[which(fiveprime_db$chr == chrom),]$start[order(abs(x - fiveprime_db[which(fiveprime_db$chr == chrom),]$start))][1])
  all.introns[which(all.introns$chr == chrom),]$end_match <- sapply(all.introns[which(all.introns$chr == chrom),]$end, function(x) threeprime_db[which(threeprime_db$chr == chrom),]$start[order(abs(x - threeprime_db[which(threeprime_db$chr == chrom),]$start))][1])

}

all.introns$start_diff <- all.introns$start - all.introns$start_match
all.introns$end_diff <- all.introns$end - all.introns$end_match

##Strand adjust the 5p and 3p distances:
all.introns$fivep_diff <- NA 
all.introns$threep_diff <- NA

all.introns[which(all.introns$strand_1 == "+"), ]$fivep_diff <- all.introns[which(all.introns$strand_1 == "+"), ]$start_diff
all.introns[which(all.introns$strand_1 == "+"), ]$threep_diff <- all.introns[which(all.introns$strand_1 == "+"), ]$end_diff
  
all.introns[which(all.introns$strand_1 == "-"),]$fivep_diff <- (all.introns[which(all.introns$strand_1 == "-"), ]$end_diff)*-1
all.introns[which(all.introns$strand_1 == "-"),]$threep_diff <- (all.introns[which(all.introns$strand_1 == "-"), ]$start_diff)*-1


##Create all the intersection databases
all.introns_intersect = all.junctions %>%
  left_join(intron_db, by=c("chr","start","end"))

#all.introns_intersect[which(all.introns_intersect$gene == "DYNLL1"),]

threeprime_intersect = all.junctions %>%
  select(chr, clusterID, start=end) %>%
  left_join(threeprime_db, by=c("chr","start"))

fiveprime_intersect =  all.junctions %>%
  select(chr, clusterID, start) %>%
  left_join(fiveprime_db, by=c("chr","start"))

print("Annotating junctions")

verdict.list <- list()
strand.list <- list()
coord.list <- list()
gene.list <- list()
ensemblID.list <- list()
transcripts.list <- list()
constitutive.list <- list()

clusters <- unique( all.junctions$clusterID)

for( clu in clusters){
  # for each intron in the cluster, check for coverage of both
  # output a vector of string descriptions
  cluster <- all.junctions %>% filter( clusterID == clu )

  # first subset the intersected files to speed up later query - this uses the data.tables method
  fprimeClu <- fiveprime_intersect %>% filter( clusterID == clu )
  tprimeClu <- threeprime_intersect %>% filter( clusterID == clu )
  bothSSClu <- all.introns_intersect %>% filter( clusterID == clu )


  # for each intron in the cluster:
  # create vector of overlapping splice sites, indexed by the row of the intersect

  # five prime splice sites
  fprime=cluster %>% left_join(fprimeClu, by=c("chr","start"))

  # three prime splice sites
  tprime=cluster %>% left_join(tprimeClu, by=c("chr"="chr","end"="start"))

  # both splice sites
  bothSS=cluster %>% left_join(bothSSClu, by=c("chr","start","end"))

  # find gene and ensemblID by the most represented gene among all the splice sites - lazy
  cluster_gene <- names(sort(table(c(tprime$gene,fprime$gene)), decreasing = TRUE ))[1]

    # if no cluster gene found then leave as "."
  if( is.null(cluster_gene) ){
    cluster_gene <- "."
  }

  gene_strand <- NA
  if( cluster_gene != "." ){
    # get strand the same way - would prefer to use the strand of the junction
    strands <- c(tprime$strand, fprime$strand)
    # hope that all junctions align to the same gene on the same strand
    gene_strand <- unique( strands[ strands != "." & !is.na(strands) ])
    if( all(is.na(gene_strand)) | length(gene_strand) != 1 ){
      gene_strand <- NA
    }
  }
# do the same for EnsemblID
  cluster_ensemblIDs <- names(sort(table( c(tprime$gene_id,fprime$gene_id)), decreasing = TRUE ))
  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]

  if( length( cluster_ensemblID ) == 0 ){
    cluster_ensemblID == "."
  }

  verdict <- c()
  coord <- c()
  gene <- c()
  strand <- c()
  ensemblID <- c()
  transcripts <- list()

  
  for( intron in 1:nrow(cluster) ){

    coord[intron] <- paste(cluster[intron,]$chr,cluster[intron,]$start, cluster[intron,]$end )
    strand[intron] <- gene_strand
    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID

    fprime_intron=cluster[intron,] %>% left_join(fprime, by=c("chr","start"))
    tprime_intron=cluster[intron,] %>% left_join(tprime, by=c("chr","end"))
    bothSS_intron=cluster[intron,] %>% left_join(bothSSClu, by=c("chr","start","end"))

    # for each intron create vector of all transcripts that contain both splice sites
    transcripts[[intron]] <- intersect( tprime_intron$transcript,fprime_intron$transcript )

    verdict[intron] <- "error"

    unknown_3p=all( is.na(tprime_intron$gene) )
    unknown_5p=all( is.na(fprime_intron$gene) )

    if (is.na(gene_strand)) {
      verdict[intron] <- "unknown_strand"
    } else {
      if( all( is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene))){
        verdict[intron] <- "cryptic_unanchored"
      }
      if( (all( is.na(tprime_intron$gene )) & all( !is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
        ( all( is.na(fprime_intron$gene )) & all( !is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_threeprime"
      }
      if(
        ( all( !is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
        ( all( !is.na(fprime_intron$gene )) & all( is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_fiveprime"
      }
      if( is.na(gene_strand) & ( all( !is.na(tprime_intron$gene )) | all( !is.na(fprime_intron$gene ) ) ) ){
        verdict[intron] <- "cryptic"
      }
      if( # if both splice sites are annotated
        all( !is.na(tprime_intron$gene ) ) & all( !is.na(fprime_intron$gene ) )
      ){
        # test if the splice sites are paired in a known intron
        if( all( !is.na(bothSS_intron$gene )) ){
          verdict[intron] <- "annotated"
        }else{ # both are annotated but never in the same junction
          verdict[intron] <- "novel annotated pair"
        }
      }
    }
    
    verdict.list[[clu]] <- verdict
    coord.list[[clu]] <- coord
    gene.list[[clu]] <- gene
    strand.list[[clu]] <- strand
    ensemblID.list[[clu]] <- ensemblID
  
    # once all the transcripts for all the introns are found, go back and work out how many constitutive each junction is. Does the junction appear in every transcript?

    if( intron == nrow(cluster)){ # only on final intron
      all_transcripts <- unique( unlist( transcripts ) )
      # remove "." - non-existent transcripts
      all_transcripts <- all_transcripts[ all_transcripts != "." ]

      constitutive <- lapply( transcripts, FUN = function(x) {
        # for each intron how many transcripts is it seen in?
        x <- x[ x != "." ]
        length(x) / length( all_transcripts)

        })

      constitutive.list[[clu]] <- constitutive

      # collapse all.introns transcripts for each intron into a single string
      transcripts.list[[clu]] <- lapply(transcripts, FUN = function(x) paste( x, collapse = "+" ) )

    }

  }

}

print("Preparing results")

# match the lists together
all.introns$strand <- unlist(strand.list)[ match( paste(all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$verdict <- unlist(verdict.list)[ match( paste(all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$gene <- unlist(gene.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$ensemblID <- unlist(ensemblID.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$transcripts <- unlist( transcripts.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

all.introns$constitutive.score <-  unlist( constitutive.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

#all.introns$prediction <-  unlist( classification.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

# replace NA values/missing transcripts with "."      ??? SHould i add the start_match and end_match to this list - should not have any missing values

all.introns %<>% mutate( gene=ifelse(is.na(gene), ".", gene),
                         ensemblID=ifelse(is.na(ensemblID), ".", ensemblID),
                         transcripts=ifelse(transcripts == "", ".", transcripts),
                         constitutive.score=signif(constitutive.score, digits = 2))


print("Summary Counts for each junc type:")
table(all.introns$verdict)
print("Total number of junctions pre-filtering:")
nrow(all.introns)

print("Saving Introns info")

all_introns_meta <- all.introns

save(all_introns_meta,version=2, file = paste0(output,"_all.introns.info.Rdata") )
write.table(all_introns_meta, file= paste0(output, "_all.introns.info.txt"), quote=F, sep="\t")

print("Done")

