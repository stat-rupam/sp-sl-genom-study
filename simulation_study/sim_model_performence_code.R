# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)
K = 50
theta = 0.5
n = 5000
folder = "simulated_data_iii/"
file_path_data = paste0("~/personal/spike-slab-analysis/reusable_functions/simulation_studies/",folder)
file_name_data = paste0("sim_data_K_",K,"_theta_",(100*theta),".RData")
file_path_model = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
file_name_model = paste0("sp_sl_model_K_",K,"_theta_",(100*theta),".RData")
file_name_data = paste0(file_path_data,file_name_data)
file_name_model = paste0(file_path_model,file_name_model)
load(file_name_model)
load(file_name_data)

#actual R_p
R_p_act = sim_data[['R_p_act']]
# Create a model object for further analysis
model <- sp_sl_model

# Extract the estimated beta coefficients and pi values from the model
beta <- model$postMeans$beta
pi <- c(1, model$postMeans$pV1)
coef <- data.frame(index = 1:length(beta), beta, pi)
coef <- coef[order(coef$pi), ]

#extracing the test set
SNPS = scale(sim_data[['SNPs']])
testing_set = SNPS[(n+1):(2*n),]
y_test = sim_data[['Y']][(n+1):(2*n)]

# Prepare the predictor variables for the testing set
X_new <- as.matrix(testing_set)
y_new <- y_test

# Initialize an empty vector for storing results
ss <- NULL
num_of_var <- NULL

# Loop through different probability thresholds
for (prob in coef$pi[-length(coef$pi)]) {
  # Filter coefficients based on the current probability threshold
  coef_new <- coef[coef$pi > prob,]
  
  # Select the corresponding predictors
  df1_test <- X_new[, coef_new$index]
  
  # Calculate predicted values using selected coefficients
  p_y_hat <- as.matrix(df1_test) %*% coef_new$beta
  
  
  # Calculate R-squared for the current model
  ss_y_hat <- 1 - (sum((y_new - p_y_hat)^2) / n) / var(y_new)
  size <- dim(coef_new)[1]
  print(ss_y_hat)
  
  # Append the R-squared value to the results vector
  ss <- c(ss, ss_y_hat)
  num_of_var <- c(size, num_of_var)
}

# Create a data frame to store R-squared values and associated probabilities
R_p <- ss
print(data.frame(R_p_act,max(R_p),K-num_of_var[R_p==max(R_p)]))
print(sp_sl_model$hyperparameters$w)
prob <- coef$pi[-length(coef$pi)]
df <- data.frame(R_p, prob, num_of_var)

