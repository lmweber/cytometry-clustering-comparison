#########################################################################################
# R script to run flowMeans
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(flowMeans)


# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))

data_flowMeans <- data[, -grep("label", colnames(data))]

head(data_flowMeans)
dim(data_flowMeans)

# run flowMeans

system.time(
res_flowMeans <- flowMeans(data_flowMeans, Standardize = FALSE)
)

save(res_flowMeans, file = "../results/Levine_BMMC_32/flowMeans/res_flowMeans.RData")

# save cluster labels

clus_flowMeans <- res_flowMeans@Label

res <- data.frame(label = clus_flowMeans)
write.table(res, file = "../results/Levine_BMMC_32/flowMeans/Levine_BMMC_32_flowMeans_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

