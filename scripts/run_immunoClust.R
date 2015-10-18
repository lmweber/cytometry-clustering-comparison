# R script to run immunoClust
# Lukas Weber, Oct 2015

# note installation from Bioconductor requires GNU Scientific Library


library(flowCore)
library(immunoClust)


# load data (non-transformed)

file <- "../data/Levine_BMMC_32_H1/Levine_BMMC_32_H1_notransform.fcs"
data <- flowCore::read.FCS(file, transformation = FALSE)

data
dim(data)

# truth labels

clus_truth <- flowCore::exprs(data)[, "label"]
table(clus_truth)
length(clus_truth)

# run immunoClust

pars <- colnames(data)[-grep("label", colnames(data))]
pars

out_immunoClust <- immunoClust::cell.process(data, parameters = pars)  # can use bias argument to tune no. of clusters

save(out_immunoClust, file = "../results/Levine_BMMC_32_H1/res_immunoClust.RData")

summary(out_immunoClust)

# plot results

data_transf <- immunoClust::trans.ApplyToData(out_immunoClust, data)

immunoClust::splom(out_immunoClust, data_transf, N = 1000)

