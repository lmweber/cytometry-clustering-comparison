#########################################################################################
# R script to prepare benchmark data set Nilsson_rare
#
# This is a 13-dimensional flow cytometry data set containing several cell populations, 
# including a rare population of hematopoietic stem cells (HSCs).
#
# This R script pre-processes the data set, adds manually gated cell population labels, 
# and exports it in .txt and .fcs formats. Gating was previously performed in Cytobank.
#
# Source: Figure 2 in the following paper:
# Nilsson et al. (2013), "Frequency Determination of Rare Populations by Flow Cytometry: 
# A Hematopoietic Stem Cell Perspective", Cytometry Part A, 83A, 721-727.
#
# Link to paper: http://www.ncbi.nlm.nih.gov/pubmed/23839904
# Link to data: http://flowrepository.org/id/FR-FCM-ZZ6L
#
# Lukas Weber, July 2016
#########################################################################################


# load packages

library(flowCore)  # from Bioconductor
library(magrittr)  # from CRAN




#################
### LOAD DATA ###
#################


# -------------
# load raw data
# -------------

# data from Figure 2 in Nilsson et al. (2013) paper (see above for details and links)

file_raw <- list.files("raw_data", pattern = "\\.fcs$", full.names = TRUE)

data_raw <- flowCore::read.FCS(file_raw, transformation = FALSE, truncate_max_range = FALSE)

head(data_raw)
dim(data_raw)  # 119,959 events, 13 protein parameters (plus scatter, live/dead, time)

parameters(data_raw)@data  # CD protein names (note PI is for dead cells)



# --------------------------------------------------------
# load data file containing single cells, non-debris, live
# --------------------------------------------------------

# single cells, non-debris, live

# pre-gating to exclude doublets, debris, and dead cells was done in Cytobank, following
# the gating scheme shown in the first 3 panels of Figure 2 in Nilsson et al. (2013)

# note compensation is done automatically in Cytobank using the spillover matrix from the
# FCS file, which can be accessed with:
# description(data_raw)$SPILL

file_single_nondeb_live <- list.files("gated_data_from_Cytobank", pattern = "Live\ cells\\.fcs$", full.names = TRUE)

data_single_nondeb_live <- flowCore::exprs(flowCore::read.FCS(file_single_nondeb_live, 
                                                              transformation = FALSE, 
                                                              truncate_max_range = FALSE))

head(data_single_nondeb_live)
dim(data_single_nondeb_live)  # 44,140 cells


# column names (protein markers)

col_names <- parameters(flowCore::read.FCS(file_single_nondeb_live))$desc
length(col_names)

col_names %>% 
  subset(., grepl("^CD", .)) %>% 
  gsub("^(CD[0-9A-Za-z]+).*", "\\1", .) %>% 
  unname -> 
  cols_markers

cols_markers
length(cols_markers)

colnames(data_single_nondeb_live) <- col_names
ix_markers <- grep("^CD", col_names)
colnames(data_single_nondeb_live)[ix_markers] <- cols_markers



# ------------------------------------------------
# load data file containing HSCs (rare population)
# ------------------------------------------------

# rare population of hematopoietic stem cells (HSCs)

# gating was done in Cytobank, following the gating scheme in Figure 2 in Nilsson et al. (2013)

file_HSC <- list.files("gated_data_from_Cytobank", pattern = "HSCs\\.fcs$", full.names = TRUE)

data_HSC <- flowCore::exprs(flowCore::read.FCS(file_HSC, transformation = FALSE, truncate_max_range = FALSE))

head(data_HSC)
dim(data_HSC)  # 358 cells

nrow(data_HSC) / nrow(data_single_nondeb_live) * 100  # 0.8% of single, non-debris, live cells




################################
### LABEL POPULATION OF HSCS ###
################################

# since the FCS files do not include event numbers or identifiers, we identify cells from
# the HSC population by checking for duplicate rows

data_dup <- rbind(data_single_nondeb_live, data_HSC)
dim(data_dup)

ix_dup <- duplicated(data_dup, fromLast = TRUE)

sum(ix_dup)  # 358 cells
length(ix_dup)

# create vector of labels

labels <- ix_dup[1:nrow(data_single_nondeb_live)]

head(labels)
sum(labels)     # 358 HSCs
length(labels)  # 44,140 single, non-debris, live cells

# add labels to data frame

data <- cbind(data_single_nondeb_live, label = labels)

head(data)
dim(data)
table(data[, "label"])




#########################
### ARCSINH TRANSFORM ###
#########################

# arcsinh transform
# using scale factor 150 for flow cytometry data (see Bendall et al. 2011, Supp. Fig. S2)

data_notransform <- data

asinh_scale <- 150

cols_to_scale <- 5:18  # columns to transform
data[, cols_to_scale] <- asinh(data[, cols_to_scale] / asinh_scale)

summary(data[, cols_to_scale])




###################
### EXPORT DATA ###
###################

# save data files in TXT and FCS format

write.table(data, file = "data/Nilsson_rare.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_notransform, file = "data/Nilsson_rare_notransform.txt", quote = FALSE, sep = "\t", row.names = FALSE)

flowCore::write.FCS(flowCore::flowFrame(data), filename = "data/Nilsson_rare.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_notransform), filename = "data/Nilsson_rare_notransform.fcs")


