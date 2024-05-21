install.packages("cvTools","BoomSpikeSlab")
#load required model configs
source("~/personal/spike-slab-analysis/reusable_functions/model_configs.R")
library(BoomSpikeSlab)
library(cvTools)

# Define the cross-validation function
cv_calculation <- function(data, response, n_fold, niter = 20000, 
                           prior_information_weight = 0.01, 
                           diagonal_shrinkage = 0.01, 
                           burn = 2000) {
  n <- dim(data)[1]
  
  # Combine data and response into a single data frame
  dataset <- data.frame(data, response)
  colnames(dataset)[ncol(dataset)] <- "y" # Rename response variable
  
  # Create k folds
  folds <- cvFolds(n = nrow(dataset), K = n_fold)
  
  # Initialize lists to store performance metrics and model summaries
  model_list <- list()
  performance <- list()
  
  # Perform k-fold cross-validation
  for (i in 1:n_fold) {
    # Split the data into training and test sets
    print(paste0("Running MCMC For The Validation Set: ",i," ,Chain Length: ", niter))
    test_indices <- which(folds$which == i)
    train_indices <- setdiff(1:nrow(dataset), test_indices)
    
    training_set <- dataset[train_indices, ]
    test_set <- dataset[test_indices, ]
    
    # Fit the model on the training set
    sp_sl_model_new <- lm.spike(
      formula = y ~ . - 1,
      data = training_set,
      niter = niter,
      prior.information.weight = prior_information_weight,
      diagonal.shrinkage = diagonal_shrinkage
    )
    
    # Store the model summary
    model_list[[i]] <- summary(object = sp_sl_model, burn = burn)$coefficients
    
    # Predict on the test set
    predictions <- predict(object = sp_sl_model, newdata = test_set, burn = burn, mean.only = TRUE)
    
    # Evaluate performance (e.g., mean squared error)
    mse <- mean((test_set$y - predictions)^2)
    
    # Store performance metrics
    performance[[i]] <- mse
  }
  
  # Summarize performance
  avg_performance <- mean(unlist(performance))
  
  # Return both model summaries and performance metrics
  return(list(models = model_list, performance = performance, avg_performance = avg_performance))
}

# Example usage
n <- 5000
training_set <- SNPS[1:n,]
test_set <- SNPS[(n+1):(2*n),]
Y_train <- Y[1:n]
Y_test <- Y[(n+1):(2*n)]

# Extract the predictors (features) and choose a specific phenotype for analysis
X <- training_set  # Extract the predictor variables
y <- Y_train # Extract the response variable

# Perform k-fold cross-validation
n_fold <- 5 # Number of folds
cv_results <- cv_calculation(data = X, response = y, n_fold = n_fold)

# Print the average performance
print(cv_results$avg_performance)
