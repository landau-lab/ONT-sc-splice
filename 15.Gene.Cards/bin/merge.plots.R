#!/usr/bin/env Rscript
#options(echo=TRUE)
## Merge Plots

library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(optparse)

option_parser=OptionParser(
  usage="%prog [options] dPSI.Exp.plot sashimi.plot output.file"
)

parsed_args <- parse_args(option_parser,  positional_arguments = 3)

dPSI.Exp.plot <- parsed_args$args[1]
sashimi.plot <- parsed_args$args[2]
output.file <- parsed_args$args[3]

img.1 <- readPNG(dPSI.Exp.plot)
g.1 <- rasterGrob(img.1, interpolate=TRUE)

plot.1 <- qplot(1:10, 1:10, geom="blank") + theme_transparent() + theme(plot.margin = margin(0.5, 0, 0, 1, "cm")) + annotation_custom(g.1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

img.2 <- readPNG(sashimi.plot)
g.2 <- rasterGrob(img.2, interpolate=TRUE)

plot.2 <- qplot(1:10, 1:10, geom="blank") + theme_transparent() + annotation_custom(g.2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

pdf(output.file, width = 100, height = 100)
plot_grid(plot.1, plot.2)
dev.off()
