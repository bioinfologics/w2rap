# Script:  nextclip_plot_duplication.R
# Purpose: R script for plotting PCR duplication distributions for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

args <- commandArgs(trailingOnly = TRUE)

input_filename <- args[1]
output_filename <- args[2]

pdf(output_filename, width=3, height=2, pointsize=6)
par(mar = c(3, 3, 1.5, 0.7))
par(mgp = c(1.8, 0.5, 0))
input_file = read.table(input_filename, nrow=100, header=TRUE)
plot(input_file[,1], input_file[,3], type="o", pch=20, col="red", ann=FALSE, cex=0.75)
title(xlab="Duplication level")
title(ylab="% read pairs")
garbage <- dev.off()

