# Script:  nextclip_plot_lengths.R
# Purpose: R script for plotting read lengths for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

args <- commandArgs(trailingOnly = TRUE)

input_filename <- args[1]
output_filename <- args[2]

pdf(output_filename, width=3, height=2, pointsize=6)
par(mar = c(3, 3, 0.27, 0.3))
par(mgp = c(1.8, 0.5, 0))
input_file = read.table(input_filename, header=FALSE)
barplot(input_file[,3], names=input_file[,1], border=NA, col="red", xlab="Read length", ylab="Number of pairs >= length", space=0)
garbage <- dev.off()

