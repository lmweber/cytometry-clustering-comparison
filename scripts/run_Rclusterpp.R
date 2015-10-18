# R script to run Rclusterpp
# Lukas Weber, Oct 2015


library(Rclusterpp)


# load data (already transformed)

file <- "../data/Levine_BMMC_32_H1/Levine_BMMC_32_H1.txt"
data <- read.table(file, header = TRUE, sep = "\t")
head(data)
dim(data)

clus_truth <- data[, "label"]
table(clus_truth)
max(clus_truth)

data_Rclusterpp <- data[, -grep("label", colnames(data))]
head(data_Rclusterpp)
dim(data_Rclusterpp)


# run Rclusterpp

out_Rclusterpp <- Rclusterpp.hclust(data_Rclusterpp, method = "complete", distance = "euclidean")

save(out_Rclusterpp, file = "../results/Levine_BMMC_32_H1/res_Rclusterpp.RData")
