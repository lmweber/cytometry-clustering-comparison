#########################################################################################
# R script to run immunoClust
#
# Lukas M. Weber, October 2015
#########################################################################################


# note installation from Bioconductor requires GNU Scientific Library


library(flowCore)
library(immunoClust)


# use non-transformed input data, since immunoClust will transform automatically

# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32_notransf.fcs"
data <- flowCore::read.FCS(file, transformation = FALSE)

data
dim(data)

# select parameters

pars <- colnames(data)[-grep("label", colnames(data))]
length(pars)


# run immunoClust
# can use bias argument to tune number of clusters

system.time(
out_immunoClust <- immunoClust::cell.process(data, parameters = pars)
)

save(out_immunoClust, file = "../results/Levine_BMMC_32/immunoClust/res_immunoClust.RData")

summary(out_immunoClust)  # number of clusters


# run immmunoClust with additional step to classify all cells ("immunoClust_all")

system.time(
out_immunoClust_all <- immunoClust::cell.process(data, parameters = pars, classify.all = TRUE)
)

save(out_immunoClust_all, file = "../results/Levine_BMMC_32/immunoClust_all/res_immunoClust_all.RData")

summary(out_immunoClust_all) # number of clusters


# save cluster assignments

clus_immunoClust <- out_immunoClust@label
res <- data.frame(label = clus_immunoClust)
write.table(res, file = "../results/Levine_BMMC_32/immunoClust/Levine_BMMC_32_immunoClust_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

clus_immunoClust_all <- out_immunoClust_all@label
res_all <- data.frame(label = clus_immunoClust_all)
write.table(res_all, file = "../results/Levine_BMMC_32/immunoClust_all/Levine_BMMC_32_immunoClust_all_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# plots

data_transf <- immunoClust::trans.ApplyToData(out_immunoClust, data)

png("../results/Levine_BMMC_32/immunoClust/Levine_BMMC_32_immunoClust.png", 
    width = 1000, height = 1000)
immunoClust::splom(out_immunoClust, data_transf, N = 1000)
dev.off()


data_transf_all <- immunoClust::trans.ApplyToData(out_immunoClust_all, data)

png("../results/Levine_BMMC_32/immunoClust_all/Levine_BMMC_32_immunoClust_all.png", 
    width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all, data_transf_all, N = 1000)
dev.off()

