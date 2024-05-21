#code for model training and model performance check
# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)
# Load necessary libraries
library(BoomSpikeSlab)  # Load the boomSpikeSlab package for the analysis
library(coda)          # Load the coda package for MCMC chain diagnostics
library(Cairo)         # Load the Cairo package for generating graphics

# Load training and testing data (temporary, for re-usability; may be removed later)
# The paths to data files are specified here
K = 5000
theta = 0.8
file = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
file_name = paste0("sim_data_K_",K,"_theta_",(100*theta),".RData")
folder = "simulated_data_vi/"
file = paste0(file,folder,file_name)
load(file)  # Load the training data

#load required functions
source("~/personal/spike-slab-analysis/reusable_functions/cv_code.R")
SNPS = scale(sim_data[['SNPs']])
Y = sim_data[['Y']]
R_p_act = sim_data[['R_p_act']]
#loading the training data
n = 5000
training_set = SNPS[1:n,]
test_set = SNPS[(n+1):(2*n),]
Y_train = Y[1:n]
Y_test = Y[(n+1):(2*n)]
# Extract the predictors (features) and choose a specific phenotype for analysis
X <- training_set  # Extract the predictor variables
y <- Y_train # Extract the response variable

cv_model <- cv_calculation(data = SNPs,
                            response = Y,
                           n_fold = 10)

# sp_sl_model <- lm.spike(formula = y~.-1,
#                         data = data.frame(X,y),
#                         niter = 300,
#                         prior.information.weight = 10^-3,
#                         diagonal.shrinkage = 10^-5)

# Get the posterior means of the model parameters
# coef <- sp_sl_model$postMeans

setwd("~/personal/spike-slab-analysis/reusable_functions/simulation_studies")
# Save the fitted model to a file for later use
file_name = paste0("sp_sl_model_K_",K,"_theta_",(100*theta),".RData")
save(sp_sl_model, file = file_name)



# Source the custom plot_function.R script
source("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis/reusable_functions/plot_function.R")


