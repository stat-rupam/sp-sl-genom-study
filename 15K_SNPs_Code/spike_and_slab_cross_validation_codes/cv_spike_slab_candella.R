install.packages("cvTools","BoomSpikeSlab")
#load required model configs
source("~/personal/spike-slab-analysis/Cluster_codes/cross_validation_codes/model_configs.R")
library(BoomSpikeSlab)
library(cvTools)

# Load necessary library
library(cvTools)
cv_analysis_single_phenotype <- function(data, phenotype, n_fold, 
                                         niter = 25000, burn = 5000, 
                                         prior_information_weight = 10^-3, 
                                         diagonal_shrinkage = 10^-3) {
  # Combine predictors and response into a single dataset
  response <- data[[phenotype]]
  dataset <- data.frame(data[ , -c(1:10)], response)
  colnames(dataset)[ncol(dataset)] <- "y" # Rename response variable
  
  cat("\n----------------------------------------\n")
  cat(paste("Starting cross-validation for Phenotype:", phenotype, "\n"))
  cat("----------------------------------------\n")
  
  # Noting the starting time
  start.time <- Sys.time()
  # Create k folds
  set.seed(1234)
  folds <- cvFolds(n = nrow(dataset), K = n_fold)
  n <- nrow(dataset)
  # Initialize lists to store model summaries, performance metrics, and test indices
  model_list <- list()
  performance <- list()
  test_indices_list <- list()
  model_optimum <- NULL  # To store the best model
  min_mse <- Inf         # Initialize minimum MSE to infinity
  optimum_cv_step <- NULL  # To store the CV step (fold index) of the optimum model
  
  # Perform k-fold cross-validation
  for (i in 1:n_fold) {
    cat("\n----------------------------------------\n")
    cat(paste("Fold", i, "of", n_fold, "for Phenotype:", phenotype, "\n"))
    cat("----------------------------------------\n")
    
    # Split the data into training and test sets
    cat("Splitting data into training and test sets...\n")
    test_indices <- which(folds$which == i)
    train_indices <- setdiff(1:nrow(dataset), test_indices)
    training_set <- dataset[train_indices, ]
    test_set <- dataset[test_indices, ]
    
    # Fit the Spike and Slab model
    cat("Fitting the model...\n")
    sp_sl_model <- lm.spike(
      formula = y ~ . - 1,
      data = training_set,
      niter = niter,
      prior.information.weight = prior_information_weight,
      diagonal.shrinkage = diagonal_shrinkage
    )
    
    # Generate predictions on the test set
    cat("Generating predictions...\n")
    predictions <- predict(object = sp_sl_model, newdata = test_set, burn = burn)
    
    # Calculate performance metrics
    cat("Calculating performance metrics...\n")
    mse <- mean((test_set$y - rowMeans(predictions))^2)
    mae <- mean(abs(test_set$y - rowMeans(predictions)))
    rmse <- sqrt(mse)
    r_squared <- 1 - (mse/ var(test_set$y))
    
    # Calculate credible intervals and predictive coverage
    credible_intervals <- apply(predictions, 1, quantile, probs = c(0.025, 0.975))
    coverage <- mean(test_set$y > credible_intervals[1, ] & test_set$y < credible_intervals[2, ]) * 100
    
    # Update the optimum model if the current model has a lower MSE
    if (mse < min_mse) {
      cat("Updating optimum model...\n")
      min_mse <- mse
      model_optimum <- sp_sl_model
      optimum_cv_step <- i  # Store the current fold index
    }
    
    # Store the model summary, performance metrics, and test indices
    model_list[[i]] <- summary(object = sp_sl_model, burn = burn)$coefficients
    performance[[i]] <- list(
      'RMSE' = rmse,
      'MSE' = mse,
      'MAE' = mae,
      'R^2' = r_squared,
      'Predictive Coverage (%)' = coverage
    )
    test_indices_list[[i]] <- test_indices
  }
  
  # Aggregate performance metrics
  performance_df <- do.call(rbind, lapply(performance, as.data.frame))
  #Noting the ending time
  end.time <- Sys.time()
  #Calculating total time taken by the system
  time.taken <- c(time.taken_ridge, end.time - start.time)
  print(paste("Total time taken for the cross validation: ",time.taken))
  
  # Return results
  return(list(
    models = model_list,
    performance = performance,
    test_indices = test_indices_list,
    model_optimum = model_optimum,
    optimum_cv_step = optimum_cv_step
  ))
}
