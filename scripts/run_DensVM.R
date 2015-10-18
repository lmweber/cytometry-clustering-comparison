# Script to run DensVM on Levine_BMMC_32_H1 data set
# DensVM is included in the cytofkit R/Bioconductor package
#
# Lukas Weber, October 2015

library(cytofkit)

# use main function "cytof_tsne_densvm"
# there is also a graphical version "cytof_tsne_densvm_GUI"

?cytof_tsne_densvm

#file <- "../data/Levine_BMMC_32_H1/Levine_BMMC_32/H1.fcs"
#res_dir <- "../results/Levine_BMMC_32_H1/DensVM"

#cytof_tsne_densvm(fcsFile = file, 
#                  resDir = res_dir, 
#                  fixedNum = 20000, 
#                  verbose = TRUE)
## doesn't work since also needs parameter (marker) names
## create a txt file with them

cytof_tsne_densvm_GUI()

## I think it is doing an extra logicle transformation - should skip since I have already done asinh
## doesn't seem possible - need to create a non-transformed data file!
