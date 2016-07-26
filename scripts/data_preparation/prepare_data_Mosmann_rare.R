#########################################################################################
# R script to prepare benchmark data set Mosmann_rare
#
# This is a 14-dimensional flow cytometry data set containing a rare population of live 
# activated (cytokine-producing) memory CD4 T cells, from healthy human peripheral blood 
# mononuclear cells (PBMCs) exposed to influenza antigens.
#
# This R script pre-processes the data set, adds manually gated cell population labels, 
# and exports it in .txt and .fcs formats. Gating was previously performed in Cytobank.
#
# Source: Figure 4 in the following paper:
# Mosmann et al. (2014), "SWIFT — Scalable Clustering for Automated Identification of 
# Rare Cell Populations in Large, High-Dimensional Flow Cytometry Datasets, Part 2: 
# Biological Evaluation", Cytometry Part A, 85A, 422-433.
#
# Link to paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4238823/
# Link to data: http://flowrepository.org/id/FR-FCM-ZZ8J
# (filename: "JMW034-J16OFVQX_G2 0o1 3_D07.fcs"; see Supplementary Information file 3 for
# full list of filenames)
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

# data from Figure 4 Mosmann et al. (2014) paper (see above for details and links)

file_raw <- list.files("raw_data", pattern = "\\.fcs$", full.names = TRUE)

data_raw <- flowCore::read.FCS(file_raw, transformation = FALSE, truncate_max_range = FALSE)

head(data_raw)
dim(data_raw)  # 1,600,017 events, 15 protein parameters and 6 scatter parameters


# --------------------------------------------
# load data file containing single, live cells
# --------------------------------------------

# single, live cells (i.e. excluding debris and doublets)

# pre-gating to exclude debris and doublets was done in Cytobank, following the gating
# scheme shown in the first two panels of Figure 4A in Mosmann et al. (2014)

# note compensation is done automatically in Cytobank using the spillover matrix from the
# FCS file, which can be accessed with:
# description(data_raw)$SPILL

file_single_live <- list.files("gated_data_from_Cytobank", pattern = "Live_cells\\.fcs$", full.names = TRUE)

data_single_live <- flowCore::exprs(flowCore::read.FCS(file_single_live, 
                                                       transformation = FALSE, 
                                                       truncate_max_range = FALSE))

head(data_single_live)
dim(data_single_live)  # 396,460 cells

# column names

parameters(flowCore::read.FCS(file_single_live))$desc[7:21] %>% 
  unname %>% 
  gsub(" (.*)$", "", .) %>% 
  gsub("/", "_", .) -> 
  cols_markers

cols_markers[c(5, 13)] <- c("GZB-SA", "CCL4")  # from Supplementary Info file 1
cols_markers
length(cols_markers)

colnames(data_single_live)[7:21] <- cols_markers


# ------------------------------------------------------------------
# load data file containing rare population of activated CD4 T cells
# ------------------------------------------------------------------

# rare population of live, activated CD4 T cells

# gating was done in Cytobank, following the gating scheme in Figure 4A in Mosmann et al. (2014)

file_activ <- list.files("gated_data_from_Cytobank", pattern = "IFNg_vs_TNFa\\.fcs$", full.names = TRUE)

data_activ <- flowCore::exprs(flowCore::read.FCS(file_activ, transformation = FALSE, truncate_max_range = FALSE))

head(data_activ)
dim(data_activ)  # 109 cells

nrow(data_activ) / nrow(data_single_live) * 100  # <0.03% of single, live cells




#############################################
### LABEL POPULATION OF ACTIVATED T CELLS ###
#############################################

# since the FCS files do not include event numbers or identifiers, we identify cells from
# the activated CD4 T cell population by checking for duplicate rows

data_dup <- rbind(data_single_live, data_activ)
dim(data_dup)

ix_dup <- duplicated(data_dup, fromLast = TRUE)  # takes 20 sec

sum(ix_dup)  # 109 cells
length(ix_dup)

# create vector of labels

labels <- ix_dup[1:nrow(data_single_live)]

head(labels)
sum(labels)     # 109 rare cells
length(labels)  # 396,460 single, live cells

# add labels to data frame

data <- cbind(data_single_live, label = labels)

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

cols_to_scale <- 7:21  # columns to transform
data[, cols_to_scale] <- asinh(data[, cols_to_scale] / asinh_scale)

summary(data[, cols_to_scale])




###################
### EXPORT DATA ###
###################

# save data files in TXT and FCS format

write.table(data, file = "data/Mosmann_rare.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_notransform, file = "data/Mosmann_rare_notransform.txt", quote = FALSE, sep = "\t", row.names = FALSE)

flowCore::write.FCS(flowCore::flowFrame(data), filename = "data/Mosmann_rare.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_notransform), filename = "data/Mosmann_rare_notransform.fcs")


