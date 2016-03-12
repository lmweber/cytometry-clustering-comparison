#########################################################################################
# Function to match cluster labels with "true" (manually gated) population labels, and
# calculate precision, recall, and F1 score
#
# Lukas M. Weber, March 2016
#########################################################################################


# arguments:
# - clus_algorithm: cluster labels from algorithm
# - clus_truth: true cluster labels
# (for both arguments: length = number of cells; names = cluster labels (integers))
helper_match_clusters_and_evaluate <- function(clus_algorithm, clus_truth) {
  
  # remove unassigned cells (NA's in clus_truth)
  unassigned <- is.na(clus_truth)
  clus_algorithm <- clus_algorithm[!unassigned]
  clus_truth <- clus_truth[!unassigned]
  if (length(clus_algorithm) != length(clus_truth)) warning("vector lengths are not equal")
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  pr_mat <- re_mat <- F1_mat <- matrix(NA, nrow = length(tbl_algorithm), ncol = length(tbl_truth))
  
  for (i in 1:length(tbl_algorithm)) {
    for (j in 1:length(tbl_truth)) {
      i_int <- as.integer(names(tbl_algorithm))[i]  # cluster number from algorithm
      j_int <- as.integer(names(tbl_truth))[j]  # cluster number from true labels
      
      true_positives <- sum(clus_algorithm == i_int & clus_truth == j_int, na.rm = TRUE)
      detected <- sum(clus_algorithm == i_int, na.rm = TRUE)
      truth <- sum(clus_truth == j_int, na.rm = TRUE)
      
      # calculate precision, recall, and F1 score
      
      precision_ij <- true_positives / detected
      recall_ij <- true_positives / truth
      F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
      
      if (F1_ij == "NaN") F1_ij <- 0
      
      pr_mat[i, j] <- precision_ij
      re_mat[i, j] <- recall_ij
      F1_mat[i, j] <- F1_ij
    }
  }
  
  # put back cluster labels
  
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- names(tbl_truth)
  
  # match labels using highest F1 score (note duplicates are allowed)
  
  labels_matched <- apply(F1_mat, 2, which.max)
  
  # precision, recall, F1 score, and number of cells for each matched cluster
  
  pr <- re <- F1 <- n_cells_matched <- rep(NA, ncol(F1_mat))
  names(pr) <- names(re) <- names(F1) <- names(n_cells_matched) <- names(labels_matched)
  
  for (i in 1:ncol(F1_mat)) {
    # use character names for column indices in case subsampling completely removes some true clusters
    pr[i] <- pr_mat[labels_matched[i], names(labels_matched[i])]
    re[i] <- re_mat[labels_matched[i], names(labels_matched[i])]
    F1[i] <- F1_mat[labels_matched[i], names(labels_matched[i])]
    
    n_cells_matched[i] <- sum(clus_algorithm == labels_matched[i], na.rm = TRUE)
  }
  
  return(list(pr = pr, re = re, F1 = F1, 
              labels_matched = labels_matched, n_cells_matched = n_cells_matched))
}

