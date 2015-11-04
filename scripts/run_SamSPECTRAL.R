#########################################################################################
# R script to run SamSPECTRAL
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(SamSPECTRAL)


# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))

data_SamSPECTRAL <- data[, -grep("label", colnames(data))]

head(data_SamSPECTRAL)
dim(data_SamSPECTRAL)


# run SamSPECTRAL

set.seed(123)

system.time(
out_SamSPECTRAL <- SamSPECTRAL(data_SamSPECTRAL, normal.sigma = 100, separation.factor = 1)
)

save(out_SamSPECTRAL, file = "../results/Levine_BMMC_32/SamSPECTRAL/res_SamSPECTRAL.RData")


# save cluster labels

clus_SamSPECTRAL <- out_SamSPECTRAL

table(clus_SamSPECTRAL)
length(clus_SamSPECTRAL)

res <- data.frame(label = clus_SamSPECTRAL)
write.table(res, file = "../results/Levine_BMMC_32/SamSPECTRAL/Levine_BMMC_32_SamSPECTRAL_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

