library(BoomSpikeSlab)
library(coda)
library(dplyr)

# Load training and testing data (temporary, for reusability; may be removed later)
# The paths to data files are specified here
load("~/personal/spike-slab-analysis/train_data_15000.RData")  # Load the training data
load("~/personal/spike-slab-analysis/test_data_15000.RData")   # Load the testing data
print("Loading data done.")

# Extract the predictors (features) and choose a specific phenotype for analysis
df1 <- training_set[,-c(1:10)]  # Extract the predictor variables
pheno = 'b_lab'                  # Specify the phenotype of interest
df1[[pheno]] <- training_set[[pheno]]    # Extract the response variable
print(paste0("Running the analysis for the Phenotype ", pheno))

# Model configuration
niter <- 25000 # chain length # please try another alternate of 20,000
burn <- 5000 # burn-in # please try another alternate of 5,000
prior_information_weight <- 10^-3
diagonal_shrinkage <- 10^-3

# Fit the Spike and Slab model using BoomSpikeSlab
spike_slab_model_b_lab <- lm.spike(b_lab ~ ., data = df1, niter = niter,
                                   prior.information.weight = prior_information_weight,
                                   diagonal.shrinkage = diagonal_shrinkage)
print("Model Training Done")

# Coefficient extraction from the model output
coef_df <- summary(spike_slab_model_b_lab)$coef
coef_df <- data.frame(coef_df)
print("Coefficient extraction from the model output done")

# Save model output
file = paste0("spike_slab_output_15K_", pheno, ".RData")
print("Prediction generation started")

# Prepare test data and generate predictions
df_test <- testing_set[,-c(1:10)]  # Extract the predictor variables
df_test[[pheno]] <- testing_set[[pheno]]
spike_slab_model_predict <- predict(spike_slab_model_b_lab, df_test, burn = burn)

print("Saving the model outputs")
save(spike_slab_model_b_lab, file = file)

print("Saving the prediction outputs")
file_pred = paste0("spike_slab_prediction_output_", pheno, ".RData")
save(spike_slab_model_predict, file = file_pred)