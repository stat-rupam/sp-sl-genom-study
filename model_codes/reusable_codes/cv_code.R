install.packages("cvTools","BoomSpikeSlab")
#load required model configs
source("~/personal/spike-slab-analysis/reusable_functions/model_configs.R")
library(BoomSpikeSlab)
library(cvTools)

# Load necessary library
library(cvTools)

# Define the cross-validation function
cv_calculation <- function(data, response, n_fold, niter = niter, 
                           prior_information_weight = prior_information_weight, 
                           diagonal_shrinkage = diagonal_shrinkage, 
                           burn = burn) {
  n <- dim(data)[1]
  
  # Combine data and response into a single data frame
  dataset <- data.frame(data, response)
  colnames(dataset)[ncol(dataset)] <- "y" # Rename response variable
  
  set.seed(1234)
  # Create k folds
  folds <- cvFolds(n = nrow(dataset), K = n_fold)
  
  # Initialize lists to store performance metrics and model summaries
  model_list <- list()
  performance <- list()
  
  # Perform k-fold cross-validation
  for (i in 1:n_fold) {
    # Split the data into training and test sets
    print(paste0("Running MCMC For The Validation Set: ",i))
    test_indices <- which(folds$which == i)
    train_indices <- setdiff(1:nrow(dataset), test_indices)
    
    training_set <- dataset[train_indices, ]
    test_set <- dataset[test_indices, ]
    
    # Fit the model on the training set
    sp_sl_model <- lm.spike(
      formula = y ~ . - 1,
      data = training_set,
      niter = niter,
      prior.information.weight = prior_information_weight,
      diagonal.shrinkage = diagonal_shrinkage
    )
    
    # Store the model summary
    model_list[[i]] <- summary(object = sp_sl_model, burn = burn)$coefficients
    
    # Predict on the test set
    predictions <- predict(object = sp_sl_model, newdata = test_set, burn = burn)
    
    # Calculate the 95% credible intervals for each column of predictions
    credible_intervals <- apply(predictions, 1, quantile, probs = c(0.025, 0.975))
    
    # Check if each true value is within the corresponding credible interval
    success <- sum(test_set$y > credible_intervals[1, ] & test_set$y < credible_intervals[2, ])
    predictions <- rowMeans(predictions)
    # Evaluate performance (e.g., mean squared error)
    mse <- mean((test_set$y - predictions)^2)
    mae <- mean(abs(test_set$y - predictions))
    R_2 <- 1-(mse)/var(test_set$y)
    pred_perct <- 100*(success/length(test_set$y))
    
    # Store performance metrics
    performance[[i]] <- list('RMSE' = sqrt(mse), 
                             'MSE' = mse,
                             'MAE' = mae, 
                             'R^2' = R_2,
                             'Predictive Coverage' = pred_perct)
  }
  
  # Convert list of lists to a data frame
  performance_df <- do.call(rbind, lapply(performance, as.data.frame))
  # Calculate average performance
  avg_performance <- colMeans(performance_df)
  
  # Return both model summaries and performance metrics
  return(list(models = model_list, performance = performance, avg_performance = avg_performance))
}

# Define the cross-validation function
cv_calculation_regularized <- function(data, response, n_fold) {
  n <- dim(data)[1]
  
  set.seed(1234)
  # Create k folds
  folds <- cvFolds(n = nrow(data), K = n_fold)
  
  # Initialize lists to store performance metrics and model summaries
  model_list_ridge <- list()
  performance_ridge <- list()
  model_list_lasso <- list()
  performance_lasso <- list()
  
  # Perform k-fold cross-validation
  for (i in 1:n_fold) {
    # Split the data into training and test sets
    print(paste0("Running LASSO Regression For The Validation Set: ",i))
    test_indices <- which(folds$which == i)
    train_indices <- setdiff(1:nrow(data), test_indices)
    
    training_set <- data[train_indices, ]
    y <- response[train_indices]
    test_set <- data[test_indices, ]
    y_test <- response[test_indices]
    
    # Convert data to matrix format
    X <- as.matrix(training_set)
    
    # Set lambda values for ridge and lasso models
    # lambda_values <- seq(0.01, 100, length = 100)
    
    # Perform cross-validation for lasso regression (alpha = 1)
    lasso_model <- glmnet(X, y, alpha = 1, family = "gaussian")
    predictions_lasso <- predict(lasso_model, 
                                 newx = as.matrix(test_set),
                                 type="response")
    
    # Store the model summary
    model_list_lasso[[i]] <- lasso_model
    
    
    # Evaluate performance (e.g., mean squared error)
    mse_lasso <- mean((y_test - predictions_lasso)^2)
    mae_lasso <- mean(abs(y_test - predictions_lasso))
    R_2_lasso <- 1-(mse_lasso)/var(y_test)
    
    # Store performance metrics
    performance_lasso[[i]] <-list('MSE' = mse_lasso,
                                  'RMSE' = sqrt(mse_lasso),
                                  'MAE' = mae_lasso,
                                  'R^2' = R2_lasso)
    
    print(paste0("Running Ridge Regression For The Validation Set: ",i))
    # Perform cross-validation for ridge regression (alpha = 0)
    ridge_model <- glmnet(X, y, alpha = 0, family = "gaussian")
    predictions_ridge <- predict(ridge_model, 
                                 newx = as.matrix(test_set),
                                 type="response")
    
    # Store the model summary
    model_list_ridge[[i]] <- ridge_model
    
    
    # Evaluate performance (e.g., mean squared error)
    mse_ridge <- mean((y_test - predictions_ridge)^2)
    mae_ridge <- mean(abs(y_test - predictions_ridge))
    R_2_ridge <- 1-(mse_ridge)/var(y_test)
    
    # Store performance metrics
    performance_ridge[[i]] <- list('MSE' = mse_ridge,
                                   'RMSE' = sqrt(mse_ridge),
                                   'MAE' = mae_ridge, 
                                   'R^2' = R2_ridge)
  }
  
  # Convert list of lists to a data frame
  performance_df_lasso <- do.call(rbind, lapply(performance_lasso, as.data.frame))
  performance_df_ridge <- do.call(rbind, lapply(performance_ridge, as.data.frame))
  
  # Calculate average performance
  avg_performance_lasso <- colMeans(performance_df_lasso)
  avg_performance_ridge <- colMeans(performance_df_ridge)
  
  # Return both model summaries and performance metrics
  return(list(
    models_ridge = model_list_ridge, 
    performance_ridge = performance_ridge, 
    avg_performance_ridge = avg_performance_ridge,
    models_lasso = model_list_lasso, 
    performance_lasso = performance_lasso, 
    avg_performance_lasso = avg_performance_lasso
  ))
}

