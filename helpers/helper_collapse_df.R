# helper functions to collapse list to data frame if some true clusters (rows) are 
# missing for some methods (fill empty values with zeros or NAs, and return as a data 
# frame); and remove entries for methods that are entirely missing

# collapse data frame with zeros as padding
collapse_df_zeros <- function(l) {
  # remove entries for missing methods
  missing <- sapply(l, is.null)
  l <- l[!missing]
  # total number of true clusters
  max_len <- max(sapply(l, length))
  # empty matrix
  m <- matrix(NA, nrow = max_len, ncol = length(l))
  colnames(m) <- names(l)
  rownames(m) <- 1:max_len
  # fill with values
  for (i in 1:length(l)) {
    m[, i][names(l[[i]])] <- l[[i]]
  }
  # set remaining NAs to 0
  m[which(is.na(m))] <- 0
  # return as data frame
  as.data.frame(m)
}


# collapse data frame with zeros as padding
collapse_df_NAs <- function(l) {
  # remove entries for missing methods
  missing <- sapply(l, is.null)
  l <- l[!missing]
  # total number of true clusters
  max_len <- max(sapply(l, length))
  # empty matrix
  m <- matrix(NA, nrow = max_len, ncol = length(l))
  colnames(m) <- names(l)
  rownames(m) <- 1:max_len
  # fill with values
  for (i in 1:length(l)) {
    m[, i][names(l[[i]])] <- l[[i]]
  }
  # return as data frame
  as.data.frame(m)
}


