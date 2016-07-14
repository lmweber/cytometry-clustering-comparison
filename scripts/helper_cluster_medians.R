#########################################################################################
# Function to calculate median expression of each marker and scale to min = 0, max = 1, 
# for a given cluster. Note that data should already be transformed, e.g. using standard
# asinh transform.
#
# Lukas Weber, July 2016
#########################################################################################


library(robustbase)


# arguments:
# - data: matrix of data (cells in rows, dimensions in columns)
# - labels: cluster labels
helper_cluster_medians <- function(data, labels) {
  
  data <- as.data.frame(data)
  
  # remove unassigned cells (NA's in labels)
  unassigned <- is.na(labels)
  labels <- labels[!unassigned]
  data <- data[!unassigned, ]
  if (length(labels) != nrow(data)) warning("lengths not equal")
  
  # note that data should already be transformed (e.g. asinh)
  
  split_data <- split(data, labels)
  split_data <- lapply(split_data, as.matrix)
  
  medians <- lapply(split_data, robustbase::colMedians)
  medians <- do.call(rbind, medians)
  
  # scale each column to min = 0, max = 1
  
  mins <- apply(medians, 2, min)
  maxs <- apply(medians, 2, max)
  
  medians_scaled <- scale(medians, center = mins, scale = maxs - mins)
  
  return(medians_scaled)
}

