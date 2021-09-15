#!/usr/bin/env Rscript
#options(echo=TRUE)
## generate the dPSI vs Expression plot

library(tidyverse)
library(optparse)

option_parser=OptionParser(
  usage="%prog [options] input.cryptic.junc.list full.path.output.file"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 2)

cryptic.junctions <- read.table(parsed_args$args[1])
output.file <- parsed_args$args[2]

expanded.df <- data.frame(do.call("rbind", strsplit(as.character(cryptic.junctions$V2), ":", fixed = TRUE)))
expanded.df$X2 <- as.numeric(paste(expanded.df$X2))
expanded.df$X3 <- as.numeric(paste(expanded.df$X3))
expanded.df$new.start <- NA
expanded.df$new.end <- NA

expanded.df[which(expanded.df$X4 == "+"),]$new.start <- expanded.df[which(expanded.df$X4 == "+"),]$X3 - 100
expanded.df[which(expanded.df$X4 == "+"),]$new.end <- expanded.df[which(expanded.df$X4 == "+"),]$X3 + 100

expanded.df[which(expanded.df$X4 == "-"),]$new.start <- expanded.df[which(expanded.df$X4 == "-"),]$X2 - 100
expanded.df[which(expanded.df$X4 == "-"),]$new.end <- expanded.df[which(expanded.df$X4 == "-"),]$X2 + 100

cryptic.junctions$new.junction.coords <- paste0(expanded.df$X1,":",expanded.df$new.start,"-", expanded.df$new.end)

write.table(cryptic.junctions , file=output.file, quote = F, row.names = F, col.names = F)


