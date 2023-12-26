drop_highly_correlated <- function(data, threshold) {
  # Compute the correlation matrix
  corr_matrix <- cor(data)
  
  # Find pairs of variables with correlation greater than the threshold
  high_corr_pairs <- which(abs(corr_matrix) > threshold & lower.tri(corr_matrix), arr.ind = TRUE)
  
  # Create a vector to keep track of variables to drop
  to_drop <- character(0)
  
  # Iterate over the pairs and add the second variable to the to_drop vector
  for (i in 1:nrow(high_corr_pairs)) {
    row_idx <- high_corr_pairs[i, 1]
    col_idx <- high_corr_pairs[i, 2]
    
    var1 <- colnames(data)[row_idx]
    var2 <- colnames(data)[col_idx]
    
    to_drop <- append(to_drop, var2)
  }
  
  # Drop the second variables in the highly correlated pairs
  data <- data[, !colnames(data) %in% to_drop]
  
  return(data)
}

