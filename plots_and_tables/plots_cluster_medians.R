#########################################################################################
# R script to generate plots of cluster medians (comparing median expression between
# detected clusters and manually gated populations)
#
# Lukas Weber, September 2016
#########################################################################################


library(flowCore)
library(pheatmap)
library(RColorBrewer)

# load main results (from main plots file)
load("main_results.RData")

# helper function
source("../helpers/helper_cluster_medians.R")




####################
### PREPARE DATA ###
####################

# only for data set Levine_32dim


# note: ACCENSE, ClusterX, DensVM, flowClust, flowMerge, immunoClust, SWIFT -- these
# methods need special treatment due to subsampling or non-FCS file formats

method_names <- names(res_all)
special_sub <- c("ClusterX", "DensVM", "flowClust", "flowMerge", "immunoClust")
special_format <- c("ACCENSE", "SWIFT")

# number of subsampled points
n_sub <- list(
  ClusterX = 100000, 
  DensVM = 100000, 
  flowClust = 10000, 
  flowMerge = 10000, 
  immunoClust = 100000
)


# load main data file (for methods with standard treatment)

file_main_Levine_32dim <- "../../../benchmark_data_sets/Levine_32dim/data/Levine_32dim.fcs"
marker_cols_Levine_32dim <- 5:36
data_main_Levine_32dim <- flowCore::exprs(flowCore::read.FCS(file_main_Levine_32dim, transformation = FALSE, truncate_max_range = FALSE))
data_main_Levine_32dim <- data_main_Levine_32dim[, marker_cols_Levine_32dim]


# load data for each method

data <- vector("list", length(method_names))
names(data) <- method_names

for (i in 1:length(data)) {
  if (!(names(data)[i] %in% c(special_sub, special_format))) {
    data[[i]] <- data_main_Levine_32dim
    
  } else if (names(data)[i] %in% special_sub) {
    # methods with subsampling: re-generate original subsampled data, using same random
    # seed from run scripts
    set.seed(123)
    ix <- sample(1:nrow(data_main_Levine_32dim), n_sub[[names(data)[i]]])
    data[[i]] <- data_main_Levine_32dim[ix, ]
    
  } else if (names(data)[i] == "ACCENSE") {
    # ACCENSE: load data directly from output file
    output_file <- "../../results_auto/ACCENSE/accense_output_Levine_32dim.csv"
    data[[i]] <- read.csv(output_file, stringsAsFactors = FALSE)
    # different column indices
    marker_cols_ACCENSE <- 6:37
    data[[i]] <- data[[i]][, marker_cols_ACCENSE]
    # require as matrix
    data[[i]] <- as.matrix(data[[i]])
    
  } else if (names(data)[i] == "SWIFT") {
    # SWIFT: load data directly from output file
    output_file <- "../../results_auto/SWIFT/Levine_32dim_notransform_subsampled.fcs"
    data[[i]] <- flowCore::exprs(flowCore::read.FCS(output_file, transformation = FALSE, truncate_max_range = FALSE))
    data[[i]] <- data[[i]][, marker_cols_Levine_32dim]
  }
}

# above is data for Levine_32dim only
data_Levine_32dim <- data

# remove missing methods for Levine_32dim
ix_remove_Levine_32dim <- names(data_Levine_32dim) %in% c("flowClust", "flowMerge", "SPADE")
data_Levine_32dim <- data[!ix_remove_Levine_32dim]


# retrieve cluster labels (previously loaded)

clus_methods <- list(
  ACCENSE = clus_ACCENSE, 
  ClusterX = clus_ClusterX, 
  DensVM = clus_DensVM, 
  FLOCK = clus_FLOCK, 
  flowClust = clus_flowClust, 
  flowMeans = clus_flowMeans, 
  flowMerge = clus_flowMerge, 
  flowPeaks = clus_flowPeaks, 
  FlowSOM = clus_FlowSOM, 
  FlowSOM_pre = clus_FlowSOM_pre, 
  immunoClust = clus_immunoClust, 
  kmeans = clus_kmeans, 
  PhenoGraph = clus_PhenoGraph, 
  Rclusterpp = clus_Rclusterpp, 
  SamSPECTRAL = clus_SamSPECTRAL, 
  SPADE = clus_SPADE, 
  SWIFT = clus_SWIFT, 
  Xshift = clus_Xshift
)

# select Levine_32dim; remove missing methods for Levine_32dim
clus_methods_Levine_32dim <- lapply(clus_methods, function(cl) cl[[1]])
clus_methods_Levine_32dim <- clus_methods_Levine_32dim[!ix_remove_Levine_32dim]


# retrieve true population labels (with subsampling where required; previously loaded)

clus_truth <- list(
  ACCENSE = clus_truth_ACCENSE, 
  ClusterX = clus_truth_ClusterX, 
  DensVM = clus_truth_DensVM, 
  FLOCK = clus_truth_FLOCK, 
  flowClust = clus_truth_flowClust, 
  flowMeans = clus_truth_flowMeans, 
  flowMerge = clus_truth_flowMerge, 
  flowPeaks = clus_truth_flowPeaks, 
  FlowSOM = clus_truth_FlowSOM, 
  FlowSOM_truth_pre = clus_truth_FlowSOM_pre, 
  immunoClust = clus_truth_immunoClust, 
  kmeans = clus_truth_kmeans, 
  PhenoGraph = clus_truth_PhenoGraph, 
  Rclusterpp = clus_truth_Rclusterpp, 
  SamSPECTRAL = clus_truth_SamSPECTRAL, 
  SPADE = clus_truth_SPADE, 
  SWIFT = clus_truth_SWIFT, 
  Xshift = clus_truth_Xshift
)

# select Levine_32dim; remove missing methods for Levine_32dim
clus_truth_Levine_32dim <- lapply(clus_truth, function(cl) cl[[1]])
clus_truth_Levine_32dim <- clus_truth_Levine_32dim[!ix_remove_Levine_32dim]




##################################
### HEATMAPS: TRUE POPULATIONS ###
##################################

# heatmaps showing true populations only


# retrieve data and true population labels (can select from any method without subsampling, e.g. k-means)
data_truth_Levine_32dim <- data_Levine_32dim[["kmeans"]]
clus_truth_Levine_32dim <- clus_truth_Levine_32dim[["kmeans"]]


# calculate cluster medians
# note values are already arcsinh-transformed, and each dimension will be scaled to min = 0, max = 1

medians_truth_Levine_32dim <- helper_cluster_medians(data_truth_Levine_32dim, clus_truth_Levine_32dim)
rownames(medians_truth_Levine_32dim) <- paste0("manually_gated_", rownames(medians_truth_Levine_32dim))


# plot heatmaps

filename_truth <- paste0("../../plots/Levine_32dim/cluster_medians/cluster_medians_heatmap_truth_Levine_32dim.pdf")

set.seed(123)
pheatmap(medians_truth_Levine_32dim, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_method = "average", 
         fontsize = 9, 
         filename = filename_truth, 
         width = 8, 
         height = 3.5)




###############################################
### HEATMAPS: CLUSTERS VS. TRUE POPULATIONS ###
###############################################

# heatmaps comparing detected clusters and true populations for each method


# calculate cluster medians
# note values are already arcsinh-transformed, and each dimension will be scaled to min = 0, max = 1

medians_Levine_32dim <- vector("list", length(data_Levine_32dim))
names(medians_Levine_32dim) <- names(data_Levine_32dim)

for (i in 1:length(medians_Levine_32dim)) {
  medians <- helper_cluster_medians(data_Levine_32dim[[i]], clus_methods_Levine_32dim[[i]])
  rownames(medians) <- paste0(names(clus_methods_Levine_32dim)[i], "_", rownames(medians))
  medians_Levine_32dim[[i]] <- medians
}


# plot heatmaps

plot_heights_Levine_32dim <- c(9, 9, 5.5, 8.5, 9.5, 6, 9.5, 14, 14, 9.5, 8.5, 9.5, 7, 14, 10)
fontsize_row_Levine_32dim <- c(8, 8, 8, 8, 8, 7, 8, 5, 5, 8, 8, 8, 8, 1, 8)

for (i in 1:length(medians_Levine_32dim)) {
  
  data_heatmap <- rbind(medians_truth_Levine_32dim, medians_Levine_32dim[[i]])
  
  annot_row <- data.frame(method = rep(c("manually_gated", names(medians_Levine_32dim)[i]), 
                                       times = c(nrow(medians_truth_Levine_32dim), nrow(medians_Levine_32dim[[i]]))))
  rownames(annot_row) <- rownames(data_heatmap)
  
  annot_colors <- c("red", "blue")
  names(annot_colors) <- c("manually_gated", names(medians_Levine_32dim)[i])
  annot_colors <- list(method = annot_colors)
  
  filename <- paste0("../../plots/Levine_32dim/cluster_medians/cluster_medians_heatmap_", 
                     names(medians_Levine_32dim)[i], "_Levine32dim.pdf")
  
  set.seed(123)
  pheatmap(data_heatmap, 
           color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           clustering_method = "average", 
           annotation_row = annot_row, 
           annotation_colors = annot_colors, 
           fontsize = 9, 
           fontsize_row = fontsize_row_Levine_32dim[i], 
           filename = filename, 
           width = 9.5, 
           height = plot_heights_Levine_32dim[i])
}



