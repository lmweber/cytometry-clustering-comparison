#########################################################################################
# R script to run DensVM
#
# DensVM is part of the cytofkit R/Bioconductor package
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(cytofkit)


# use main function "cytof_tsne_densvm"
# graphical version can also be launched with "cytof_tsne_densvm_GUI()"


# use non-transformed input data, since DensVM will transform automatically

# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32_notransf.fcs"
data <- flowCore::read.FCS(file, transformation = FALSE)

# note DensVM also requires a copy of the input FCS file in the output directory

out_dir <- "../results/Levine_BMMC_32/DensVM"

# select parameters

para <- colnames(data)[-grep("label", colnames(data))]
length(para)

# run DensVM
# results are saved to files in output directory

system.time(
cytofkit::cytof_tsne_densvm(fcsFile = file, resDir = out_dir, para = para, 
                            fixedNum = 20000, verbose = TRUE)
)

