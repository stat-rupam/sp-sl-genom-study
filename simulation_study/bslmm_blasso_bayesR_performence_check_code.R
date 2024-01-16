# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)
K = 1000
theta = 0.8
n = 5000
folder = "simulated_data_v/"
file_path_data = paste0("~/personal/spike-slab-analysis/reusable_functions/simulation_studies/",folder)
file_name_data = paste0("sim_data_K_",K,"_theta_",(100*theta),".RData")
file_path_model = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
model <- "bslmm" #"blasso","bslmm","bayesR"
file_name_model = paste0(model,"_model_K_",K,"_theta_",(100*theta),".RData")
file_name_data = paste0(file_path_data,file_name_data)
file_name_model = paste0(file_path_model,file_name_model)
load(file_name_model)
load(file_name_data)

#actual R_p
R_p_act = sim_data[['R_p_act']]
# Create a model object for further analysis
model <- bslmm_model

# Extract the estimated beta coefficients  from the model
beta <- model$beta

#extracing the test set
SNPS = scale(sim_data[['SNPs']])
testing_set = SNPS[(n+1):(2*n),]
y_test = sim_data[['Y']][(n+1):(2*n)]

# Prepare the predictor variables for the testing set
X_new <- as.matrix(testing_set)
y_new <- y_test

p_y_hat <- beta[1] + X_new %*% beta[2:(K+1)]


# Calculate R-squared for the current model
ss_y_hat <- 1 - (sum((y_new - p_y_hat)^2) / n) / var(y_new)

print(data.frame(ss_y_hat,R_p_act))
