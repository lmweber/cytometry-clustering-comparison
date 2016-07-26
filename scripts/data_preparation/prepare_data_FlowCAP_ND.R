#########################################################################################
# R script to prepare benchmark data set FlowCAP_ND
#
# This is a 10-dimensional flow cytometry data set from the FlowCAP-I competition 
# ("normal donors" data set). See the FlowCAP paper (reference below) for details; see
# Supplementary Table 4, and additional information at: 
# http://flowcap.flowsite.org/download/flowCAP/DataDescription.pdf
#
# This R script pre-processes the data set, adds manually gated cell population labels, 
# and exports it in .txt and .fcs formats.
#
# Note that the raw data sets are already transformed "to linear space" (see paper), so
# no additional arcsinh transformation is required.
#
# Source: FlowCAP-I competition, "normal donors" data set:
# - Aghaeepour et al. (2013), "Critical assessment of automated flow cytometry data
# analysis techniques", Nature Methods, 10(3).
# - additional information: http://flowcap.flowsite.org/
# - link to paper: http://www.nature.com/nmeth/journal/v10/n3/full/nmeth.2365.html
# - data set available on FlowRepository: https://flowrepository.org/experiments/32 or
# https://flowrepository.org/id/FR-FCM-ZZYZ (or search for "FlowCAP" in FlowRepository)
# - manual gating population labels available in Cytobank:
# https://community.cytobank.org/cytobank/experiments/4329
#
# Lukas Weber, July 2016
#########################################################################################


library(flowCore)  # from Bioconductor



#################
### LOAD DATA ###
#################

# see above for links to download files
# one FCS file per sample; separate file with manual gating population labels


# filenames

files_data <- list.files("raw_data/FlowRepository_FR-FCM-ZZYZ_files", 
                         pattern = "\\.fcs$", full.names = TRUE)

files_labels <- list.files("population_labels_from_Cytobank/FlowCAP-I/Data/Labels/NDD", 
                           pattern = "\\.csv$", full.names = TRUE)


# load FCS files; add population labels and sample numbers

data <- data.frame()

for (i in 1:length(files_data)) {
  
  data_i <- flowCore::exprs(flowCore::read.FCS(files_data[i], transformation = FALSE, truncate_max_range = FALSE))
  
  # population labels (0 = outliers, i.e. unassigned cells)
  labels_i <- read.csv(files_labels[i], header = TRUE)
  data_i <- cbind(data_i, label = labels_i[, 1])
  
  # sample number
  data_i <- cbind(data_i, sample = i)
  
  data <- rbind(data, data_i)
}

# change labels for outliers/unassigned cells to NA
data[, "label"][data[, "label"] == 0] <- NA


head(data)
tail(data)

dim(data)  # 1,778,883 cells, 10 dimensions (plus 4 other columns)
table(data[, "label"], exclude = NULL)  # 7 manual gating populations (NA = outliers/unassigned)
table(data[, "sample"])  # 30 samples


# arcsinh transform: not required, since raw data already transformed "to linear space"
# (see FlowCAP paper, Aghaeepour et al. 2013)



###################
### EXPORT DATA ###
###################

# save data file in .txt format

write.table(data, file = "data/FlowCAP_ND.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# save data file in .fcs format

flowCore::write.FCS(flowCore::flowFrame(as.matrix(data)), filename = "data/FlowCAP_ND.fcs")


