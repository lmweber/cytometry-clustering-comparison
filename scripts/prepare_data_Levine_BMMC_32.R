#########################################################################################
# R script to prepare data set "Levine_BMMC_32"
#
# Read FCS files, generate labels for manually gated clusters, apply asinh transform, and
# export in TXT and FCS format. The exported files are then ready for input to clustering
# algorithms.
#
# Source:
# Levine et al. (2015), "benchmark data set 2", 32-dimensional mass cytometry data set 
# from healthy human bone marrow samples, from two individuals H1 and H2
# - link to paper: http://www.sciencedirect.com/science/article/pii/S0092867415006376
# - link to data: https://www.cytobank.org/cytobank/experiments/46102
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(magrittr)


# read FCS filenames

# one FCS file per manually gated cluster, per individual, for two individuals H1 and H2

DIR <- "../raw_data/Levine_BMMC_32"
ds <- c("H1\\.fcs$", "H2\\.fcs$")
names(ds) <- c("H1", "H2")

files <- lapply(ds, function(p) {
  list.files(DIR, pattern = p, full.names = TRUE)
})

files

# remove "unassigned" groups (cells without cluster labels)

files <- lapply(files, function(f) {
  f <- f[-grep("NotDebrisSinglets", f)]
})

files

# check files (clusters) are in same order

names_H1 <- gsub("(^.*normalized_)|(_H1.fcs$)", "", files$H1)
names_H2 <- gsub("(^.*normalized_)|(_H2.fcs$)", "", files$H2)
names_H1
names_H2
all.equal(names_H1, names_H2)


# create column names (markers)

cols <- files$H1[1] %>% 
  read.FCS(., transformation = FALSE) %>% 
  exprs %>% 
  colnames %>% 
  unname
cols <- cols[5:36]  # remove non-protein columns

cols
length(cols)


# read FCS files individually and generate cluster labels from file number

# raw data consists of one FCS file per manually gated cluster, per individual

data <- list(matrix(nrow = 0, ncol = length(cols)), matrix(nrow = 0, ncol = length(cols)))
labels <- list(vector(), vector())
names(data) <- names(labels) <- names(files)

for (h in names(data)) {
  
  for (i in 1:length(files[[h]])) {
    
    data_i <- flowCore::exprs(flowCore::read.FCS(files[[h]][i], transformation = FALSE))
    if (!all(cols %in% colnames(data_i))) stop("column names do not match")
    data_i <- data_i[, cols]
    data[[h]] <- rbind(data[[h]], data_i)
    
    labels_i <- rep(i, nrow(data_i))
    labels[[h]] <- c(labels[[h]], labels_i)
  }
}

lapply(data, head)
lapply(data, dim)

lapply(labels, table)  # manually gated cluster labels
lapply(labels, length)


# apply asinh transformation

data_notransf <- data  # keep non-transformed data since some methods require it

asinh_scale <- 5
data <- lapply(data, function(d) {
  asinh(d / asinh_scale)
})

lapply(data, head)


# combine data and labels from individuals H1 and H2 (previously checked that clusters
# are in the same order)

res <- rbind(data$H1, data$H2)
res_notransf <- rbind(data_notransf$H1, data_notransf$H2)
res_labels <- c(labels$H1, labels$H2)

dim(res)
dim(res_notransf)
length(res_labels)


# save as tab-delimited text files

write.table(cbind(res, res_labels), 
            file = "../data/Levine_BMMC_32/Levine_BMMC_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cbind(res_notransf, res_labels), 
            file = "../data/Levine_BMMC_32/Levine_BMMC_32_notransf.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save as FCS files

flowCore::write.FCS(flowCore::flowFrame(cbind(res, res_labels)), 
                    filename = "../data/Levine_BMMC_32/Levine_BMMC_32.fcs")

flowCore::write.FCS(flowCore::flowFrame(cbind(res_notransf, res_labels)), 
                    filename = "../data/Levine_BMMC_32/Levine_BMMC_32_notransf.fcs")

