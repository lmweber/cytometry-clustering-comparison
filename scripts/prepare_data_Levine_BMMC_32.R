# R script to prepare data sets "Levine_BMMC_32_H1" and "Levine_BMMC_32_H2" from FCS 
# files. The processed files are then ready for input to clustering algorithms.
# 
# Source: Levine et al (2015), "benchmark data set 2".
# link to paper: http://www.sciencedirect.com/science/article/pii/S0092867415006376
# link to data: https://www.cytobank.org/cytobank/experiments/46102
# 
# Lukas M. Weber, October 2015


library(flowCore)
library(magrittr)

# read FCS filenames

# there is one FCS file per manually gated cluster, per individual
# and there are two individuals, H1 and H2

DIR <- "../raw_data/Levine_BMMC_32"
ds <- c("H1\\.fcs$", "H2\\.fcs$")
names(ds) <- c("H1", "H2")

files <- lapply(ds, function(p) {
  list.files(DIR, pattern = p, full.names = TRUE)
})

files

# put "unassigned" group (cells without cluster labels) at the end

files <- lapply(files, function(f) {
  f <- c(f[-grep("NotDebrisSinglets", f)], f[grep("NotDebrisSinglets", f)])
})

files

# create column names (markers)

cols <- files$H1[1] %>% 
  read.FCS(., transformation = FALSE) %>% 
  exprs %>% 
  colnames %>% 
  unname
cols <- cols[5:36]  # remove non-protein columns

cols
length(cols)

# read FCS files and generate cluster labels
# raw data consists of one FCS file per cluster, per individual

data <- list(matrix(nrow = 0, ncol = length(cols)), 
             matrix(nrow = 0, ncol = length(cols)))
labels <- list(vector(), vector())
names(data) <- names(files)
names(labels) <- names(files)

for (h in names(data)) {
  
  for (i in seq_along(files[[h]])) {
    
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
lapply(labels, table)  # these are the manually gated cluster labels for each data set
lapply(labels, length)

# apply asinh transform

asinh_scale <- 5
data <- lapply(data, function(d) {
  asinh(d / asinh_scale)
})

lapply(data, head)

# save as tab-delimited text files

res_H1 <- cbind(data[["H1"]], label = labels[["H1"]])
res_H2 <- cbind(data[["H2"]], label = labels[["H2"]])

head(res_H1)
head(res_H2)

write.table(res_H1, file = "../data/Levine_BMMC_32_H1/Levine_BMMC_32_H1.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_H2, file = "../data/Levine_BMMC_32_H2/Levine_BMMC_32_H2.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save as FCS files

res_FCS_H1 <- flowCore::flowFrame(res_H1)
res_FCS_H2 <- flowCore::flowFrame(res_H2)

flowCore::write.FCS(res_FCS_H1, filename = "../data/Levine_BMMC_32_H1/Levine_BMMC_32_H1.fcs")
flowCore::write.FCS(res_FCS_H2, filename = "../data/Levine_BMMC_32_H2/Levine_BMMC_32_H2.fcs")

