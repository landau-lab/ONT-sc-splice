#!/usr/bin/env Rscript
#options(echo=TRUE)
## Merge Plots

library(multipanelfigure)
library(optparse)

option_parser=OptionParser(
  usage="%prog [options] dPSI.Exp.plot sashimi.plot output.file"
)

dPSI.Exp.plot <- parsed_args$args[1]
sashimi.plot <- parsed_args$args[2]
output.file <- parsed_args$args[3]

figure <- multi_panel_figure(width = c(700,550), unit = "mm" , height = 1200 ,rows=1, row_spacing = 0, column_spacing =0) 

figure %<>% 
  fill_panel(dPSI.Exp.plot, row=1, column = 1) %>%
  fill_panel(sashimi.plot, row=1, column = 2) %>%
  save_multi_panel_figure(filename = output.file)

