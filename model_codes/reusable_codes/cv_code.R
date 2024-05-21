install.packages("cvTools","BoomSpikeSlab")
library(BoomSpikeSlab)
library(cvTools)

cv_calculation(data, response, n_fold){
  n <- dim(data)[1]
  # Combine data and response into a single data frame
  dataset <- data.frame(data, response)
  colnames(dataset)[ncol(dataset)] <- "y" # Rename response variable
  
  # Create k folds
  folds <- cvFolds(n = nrow(dataset), K = k)
  # Initialize a list to store performance metrics, models
  model_list <- list()
  performnce <- list()
  # Perform k-fold cross-validation
  for (i in 1:k) {
    # Split the data into training and test sets
    test_indices <- which(folds$which == i)
    train_indices <- setdiff(1:nrow(dataset), test_indices)
    
    training_set <- dataset[train_indices, ]
    test_set <- dataset[test_indices, ]
    
    # Fit the model on the training set
    sp_sl_model <- lm.spike(
      formula = y ~ . - 1,
      data = training_set,
      niter = niter,
      prior.information.weight = prior.information.weight,
      diagonal.shrinkage = diagonal.shrinkage
    )
    model_list[[i]] <- sp_sl_model$postMeans
    
    # Predict on the test set
    predictions <- predict(object = sp_sl_model, 
                           newdata = test_set,
                           burn = burn,
                           mean.only = TRUE)
    
    # Evaluate performance (e.g., mean squared error)
    mse <- mean((test_set$y - predictions)^2)
    
    # Store performance metrics
    performance[[i]] <- mse
  }
  
  # Summarize performance
  avg_performance <- mean(unlist(performance))
  return(model_list, list(performance = performance, avg_performance = avg_performance))
}