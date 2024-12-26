library(glmnet)

# Load training and testing data
load("~/personal/spike-slab-analysis/Cluster_codes/train_data_15000.RData")
load("~/personal/spike-slab-analysis/Cluster_codes/test_data_15000.RData")

# Initialize variables
phenotypes <- c()
ridge_Mse <- c()
ridge_Rmse <- c()
ridge_Mae <- c()
ridge_R2 <- c()
time.taken_ridge <- c()

lasso_Mse <- c()
lasso_Rmse <- c()
lasso_Mae <- c()
lasso_R2 <- c()
time.taken_lasso <- c()

pheno <- 'L_lab'  # Specify phenotype
training_set[[pheno]] <- as.numeric(training_set[[pheno]])
testing_set[[pheno]] <- as.numeric(testing_set[[pheno]])

# Scale predictors only
df1 <- scale(training_set[, -c(1:10)])  # Exclude metadata columns
df2 <- scale(testing_set[, -c(1:10)])

X_train <- as.matrix(df1)
y_train <- as.matrix(training_set[[pheno]])
X_test <- as.matrix(df2)
y_test <- as.matrix(testing_set[[pheno]])

# Noting the starting time
start.time <- Sys.time() 

# Perform cross-validation for Ridge regression (alpha = 0)
cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "gaussian")
best_lambda_ridge <- cv_ridge$lambda.min

# Fit Ridge model with the best lambda
ridge_model <- glmnet(X_train, y_train, alpha = 0, family = "gaussian", lambda = best_lambda_ridge)
predictions_ridge <- predict(ridge_model, newx = X_test)

# Calculate Ridge performance metrics
ridge_mse <- mean((y_test - predictions_ridge)^2)
ridge_rmse <- sqrt(ridge_mse)
ridge_mae <- mean(abs(y_test - predictions_ridge))
ridge_r2 <- 1 - sum((y_test - predictions_ridge)^2) / sum((y_test - mean(y_test))^2)

ridge_Mse <- c(ridge_Mse, ridge_mse)
ridge_Rmse <- c(ridge_Rmse, ridge_rmse)
ridge_Mae <- c(ridge_Mae, ridge_mae)
ridge_R2 <- c(ridge_R2, ridge_r2)

#Noting the ending time
end.time <- Sys.time()
#Calculating total time taken by the system
time.taken_ridge <- c(time.taken_ridge, end.time - start.time)

######------#####-----#####------######------######------#######-------########

# Noting the starting time
start.time <- Sys.time()

# Perform cross-validation for Lasso regression (alpha = 1)
cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian")
best_lambda_lasso <- cv_lasso$lambda.min

# Fit Lasso model with the best lambda
lasso_model <- glmnet(X_train, y_train, alpha = 1, family = "gaussian", lambda = best_lambda_lasso)
predictions_lasso <- predict(lasso_model, newx = X_test)

# Calculate Lasso performance metrics
lasso_mse <- mean((y_test - predictions_lasso)^2)
lasso_rmse <- sqrt(lasso_mse)
lasso_mae <- mean(abs(y_test - predictions_lasso))
lasso_r2 <- 1 - sum((y_test - predictions_lasso)^2) / sum((y_test - mean(y_test))^2)

lasso_Mse <- c(lasso_Mse, lasso_mse)
lasso_Rmse <- c(lasso_Rmse, lasso_rmse)
lasso_Mae <- c(lasso_Mae, lasso_mae)
lasso_R2 <- c(lasso_R2, lasso_r2)

phenotypes <- c(phenotypes, pheno)
#Noting the ending time
end.time <- Sys.time()
#Calculating total time taken by the system
time.taken_lasso <- c(time.taken_lasso,end.time - start.time)

# Save the models
save(ridge_model, file = paste0("ridge_model_", pheno, ".RData"))
save(lasso_model, file = paste0("lasso_model_", pheno, ".RData"))

# Store performance metrics in a data frame
performance <- data.frame(
  Phenotypes = phenotypes,
  Ridge_MSE = ridge_Mse,
  Ridge_RMSE = ridge_Rmse,
  Ridge_MAE = ridge_Mae,
  Ridge_R2 = ridge_R2,
  time_taken_ridge = time.taken_ridge,
  Lasso_MSE = lasso_Mse,
  Lasso_RMSE = lasso_Rmse,
  Lasso_MAE = lasso_Mae,
  Lasso_R2 = lasso_R2,
  time_taken_lasso = time.taken_lasso
)

# Save the performance metrics
write.csv(performance, file = paste0("performance_", pheno, ".csv"), row.names = FALSE)

# Save the models
save(ridge_model, file = paste0("ridge_model_", pheno, ".RData"))
save(lasso_model, file = paste0("lasso_model_", pheno, ".RData"))

# Display the performance table
print(performance)
