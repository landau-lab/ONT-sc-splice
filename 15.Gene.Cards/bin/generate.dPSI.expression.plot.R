#!/usr/bin/env Rscript
#options(echo=TRUE)
## generate the dPSI vs Expression plot

library(tidyverse)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(Seurat)
library(optparse)


option_parser=OptionParser(
  usage="%prog [options] gene.name cryptic.junction.coordinates path/to/bulk/alt3pss.usage.table/ path/to/ct/alt3pss.usage.tables/ celltype.list seurat.object cell.assignment.column genotype.column sample.id sample.bulk.psi.table sample.bulk.psi.column dnmt3a.bulk.psi.table"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 12)

gene <-  parsed_args$args[1]
cryptic.junction <- parsed_args$args[2] 
path.to.bulk.results <- parsed_args$args[3] 
path.to.ct.results <- parsed_args$args[4]
celltype.list <- unlist(strsplit(parsed_args$args[5],split=","))
seurat.object <- readRDS(parsed_args$args[6])
cell.assignment.column <- parsed_args$args[7]
genotype.column <- parsed_args$args[8]
sample.id <- parsed_args$args[9]
sample.bulk.psi.table <- read.table(parsed_args$args[10],sep="\t")
sample.bulk.psi.column <- parsed_args$args[11]
dnmt3a.bulk.psi.table <- read.table(parsed_args$args[12],sep="\t")

## --------------- Create dPSI vs expression results table --------------- ##

colname.list <- c("pvalue","dPSI","total.reads", "mut.psi","mut.reads", "wt.psi", "wt.reads")
psi.data <- as.data.frame(matrix(ncol=length(colname.list), nrow=1))
colnames(psi.data) <- colname.list

## --------------- Grab the results from Bulk alt.3p.ss table --------------- ##

rownames(psi.data) <- c("ALL")
cryptic.3p.splicing.table <- read.table(paste0(path.to.bulk.results,"logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt"))

psi.data["ALL",]$pvalue <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$pvalue
psi.data["ALL",]$dPSI <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$dPSI
psi.data["ALL",]$mut.psi <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$mut.psi
psi.data["ALL",]$wt.psi <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$wt.psi
psi.data["ALL",]$mut.reads <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$three.mut.cluster.cov
psi.data["ALL",]$wt.reads <- cryptic.3p.splicing.table[which(cryptic.3p.splicing.table$intron_junction == cryptic.junction),]$three.wt.cluster.cov
psi.data["ALL",]$total.reads <- psi.data["ALL",]$mut.reads + psi.data["ALL",]$wt.reads

## --------------- Grab the results from CT alt.3p.ss tables --------------- ##

for (ct in celltype.list){
  
  cryptic.3p.splicing.table.ct <- read.table(paste0(path.to.ct.results, ct ,"/logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt"))
  
  if(nrow(cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),])>0){
    psi.data[ct,]$pvalue <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$pvalue
    psi.data[ct,]$dPSI <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$dPSI
    psi.data[ct,]$mut.psi <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$mut.psi
    psi.data[ct,]$wt.psi <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$wt.psi
    psi.data[ct,]$mut.reads <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$three.mut.cluster.cov
    psi.data[ct,]$wt.reads <- cryptic.3p.splicing.table.ct[which(cryptic.3p.splicing.table.ct$intron_junction == cryptic.junction),]$three.wt.cluster.cov
    psi.data[ct,]$total.reads <- psi.data[ct,]$mut.reads + psi.data[ct,]$wt.reads
    
  }
}


celltype.list.full <- c("ALL", celltype.list)
psi.data <- add_column(psi.data, rownames(psi.data), .before = 1)
colnames(psi.data)[1] <- "CellType"
psi.data$CellType <- factor(psi.data$CellType, levels = celltype.list.full)


## --------------- Grab Normalized Expression Values from Seurat object - ONT --------------- ##


Idents(seurat.object) <- cell.assignment.column
expression.values.ct <- AverageExpression(seurat.object, assay = "RNA" ,features = gene, slot="data")
expression.values.ct$RNA <- as.data.frame(t(expression.values.ct$RNA))

seurat.object@meta.data$All <- "ALL"
Idents(seurat.object) <- "All"
expression.values.all <- AverageExpression(seurat.object, assay = "RNA" ,features = gene, slot="data")
expression.values.all$RNA <- as.data.frame(t(expression.values.all$RNA))

expression.values <- rbind(expression.values.all$RNA, expression.values.ct$RNA)

psi.data$ONT.Expression <- NA
psi.data$ONT.Expression <- expression.values[rownames(psi.data),gene]


## --------------- Generate dPSI vs Expression Plot  --------------- ##

OldMax = max(psi.data$ONT.Expression)
OldMin  = min(psi.data$ONT.Expression)
OldRange = (OldMax - OldMin)
NewMax = max(psi.data$dPSI)
NewMin = min(psi.data$dPSI)
NewRange = (NewMax - NewMin)

ExpressionColor <- "red"
dPSIColor <- "blue"

plot <- ggplot(psi.data, aes(x=CellType, group = 1)) + 
  geom_line( aes(y=dPSI), color=dPSIColor) + 
  geom_line( aes(y=(((ONT.Expression - OldMin) * NewRange) / OldRange) + NewMin) , linetype = "dashed", color=ExpressionColor)+
  scale_y_continuous(
    name = "dPSI",
    sec.axis = sec_axis(~(((. - NewMin ) * OldRange ) / NewRange ) + OldMin, name="ONT Expression \n (log-normalized counts)")) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.y = element_text(color =dPSIColor , size=13), 
        axis.title.y.right = element_text(color = ExpressionColor, size=13)) +
  ggtitle(paste0("dPSI and Expression for ", gene ," cryptic junction - ", cryptic.junction)) 


## --------------- Grab bulk PSI info  --------------- ##

sample.bulk.psi.val <- round(sample.bulk.psi.table[cryptic.junction,sample.bulk.psi.column],3)
rownames(dnmt3a.bulk.psi.table) <- dnmt3a.bulk.psi.table$intron.junction
dnmt3a.bulk.psi.val <- round(dnmt3a.bulk.psi.table[cryptic.junction,]$bulk.psi,3)

print(sample.bulk.psi.val)
print(dnmt3a.bulk.psi.val)

## --------------- Generate final plot with all info  --------------- ##
text1 <- ggparagraph(text = paste0("BULK_WT psi = ", dnmt3a.bulk.psi.val,"\n \nBULK_MUT psi = ",sample.bulk.psi.val), size = 15, color = "black") 
text1 <- plot_grid(NULL,NULL,text1, nrow = 2, ncol=3, rel_widths = c(2, 2, 1.5), rel_heights = c(2, 1))

text2 <- ggparagraph(text = paste0(sample.id, " UMI \n",gene,"\n",cryptic.junction), face = "bold", size = 20, color = "black")

table <- ggtexttable(psi.data,row = NULL,theme = ttheme("light"))


final.plot <- ggarrange(text1, text2, plot, table,NULL,
          ncol = 1, nrow = 5,
          heights = c(0.7,0.3,1,0.5,0.7))

png(paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/15.Gene.Cards/dPSI.Expression.plots/",sample.id,".",gene,".",cryptic.junction,".dPSI.Exp.plot.png"), width = 250, height = 450, units = 'mm', res = 300)
final.plot
dev.off()
