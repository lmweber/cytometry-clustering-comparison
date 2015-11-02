#########################################################################################
# R script to run Rclusterpp
#
# Lukas M. Weber, October 2015
#########################################################################################


# too slow for laptop - paste code into R session on multi-core server


library(flowCore)
library(Rclusterpp)


# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))

data_Rclusterpp <- data[, -grep("label", colnames(data))]

head(data_Rclusterpp)
dim(data_Rclusterpp)


# run Rclusterpp
# note: uses maximum number of cores if setThreads is left as default

Rclusterpp.setThreads(8)  # number of cores

system.time(
out_Rclusterpp <- Rclusterpp.hclust(data_Rclusterpp, method = "average", distance = "euclidean")
)

save(out_Rclusterpp, file = "../results/Levine_BMMC_32/Rclusterpp/res_Rclusterpp.RData")


# cut dendrogram at an arbitrary number of clusters k and extract labels

k <- 20

clus_Rclusterpp <- cutree(out_Rclusterpp, k = k)

table(clus_Rclusterpp)

res <- data.frame(label = clus_Rclusterpp)
write.table(res, file = "../results/Levine_BMMC_32/Rclusterpp/Levine_BMMC_32_Rclusterpp_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

