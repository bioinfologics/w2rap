# Script:  nextclip_plot_inserts.R
# Purpose: R script for plotting insert size distributions for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

args <- commandArgs(trailingOnly = TRUE)

input_filename <- args[1]
output_filename <- args[2]

pdf(output_filename, width=3, height=2, pointsize=6)
par(mar = c(3, 3, 1.5, 0.7))
par(mgp = c(1.8, 0.5, 0))
input_file = read.table(input_filename, header=FALSE)
hist(input_file[,2], breaks=500, freq=TRUE, border=NA, col="red", xlab="Distance between reads", ylab="Number of pairs", main=NULL)
garbage <- dev.off()

