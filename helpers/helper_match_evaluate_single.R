#########################################################################################
# Function to match cluster labels with manual gating (reference standard) population 
# labels and calculate precision, recall, and F1 score
#
# Matching criterion: maximum F1 score
#
# Use this function for data sets with a single (e.g. rare) population of interest
#
# Lukas Weber, August 2016
#########################################################################################


# arguments:
# - clus_algorithm: cluster labels from algorithm
# - clus_truth: true cluster labels (1 = rare cluster of interest, 0 = all others)
# (for both arguments: length = number of cells; names = cluster labels (integers))
helper_match_evaluate_single <- function(clus_algorithm, clus_truth) {
  
  # number of detected clusters
  n_clus <- length(table(clus_algorithm))
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  pr_mat <- re_mat <- F1_mat <- matrix(NA, nrow = length(tbl_algorithm), ncol = 1)
  
  for (i in 1:length(tbl_algorithm)) {
    i_int <- as.integer(names(tbl_algorithm))[i]  # cluster number from algorithm
    
    j_int <- 1  # true cluster number of the rare population of interest
    
    true_positives <- sum(clus_algorithm == i_int & clus_truth == j_int, na.rm = TRUE)
    detected <- sum(clus_algorithm == i_int, na.rm = TRUE)
    truth <- sum(clus_truth == j_int, na.rm = TRUE)
    
    # calculate precision, recall, and F1 score
    precision_ij <- true_positives / detected
    recall_ij <- true_positives / truth
    F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
    
    if (F1_ij == "NaN") F1_ij <- 0
    
    pr_mat[i, j_int] <- precision_ij
    re_mat[i, j_int] <- recall_ij
    F1_mat[i, j_int] <- F1_ij
  }
  
  # put back cluster labels (note some row names may be missing due to removal of unassigned cells)
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- "1"  # one column only
  
  # match label (single cluster only) using highest F1 score
  # use row names since some labels may have been removed due to unassigned cells
  labels_matched <- as.numeric(rownames(F1_mat)[apply(F1_mat, 2, which.max)])
  names(labels_matched) <- "1"  # one column only
  
  # precision, recall, F1 score, and number of cells for single matched cluster
  # use character names for row and column indices in case subsampling completely removes some clusters
  pr <- pr_mat[as.character(labels_matched), "1"]
  re <- re_mat[as.character(labels_matched), "1"]
  F1 <- F1_mat[as.character(labels_matched), "1"]
  
  n_cells_matched <- sum(clus_algorithm == labels_matched, na.rm = TRUE)
  
  return(list(n_clus = n_clus, 
              pr = pr, 
              re = re, 
              F1 = F1, 
              labels_matched = labels_matched, 
              n_cells_matched = n_cells_matched))
}


