#code for model training and model performence check
# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)
# Load necessary libraries
library(spikeSlabGAM)  # Load the spikeSlabGAM package for the analysis
library(coda)          # Load the coda package for MCMC chain diagnostics
library(Cairo)         # Load the Cairo package for generating graphics

# Load training and testing data (temporary, for reusability; may be removed later)
# The paths to data files are specified here
K = 1000
theta = 0.8
file = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
file_name = paste0("sp_sl_model_",K,"_theta_",(100*theta),".RData")
file = paste0(file,file_name)
load(file)  # Load the training data
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

# Define hyperparameters for the spike-and-slab regression model
# Hyperparameters for the Beta-Beta prior for weights (w)
# Default values are for a uniform distribution.
w <- c(alphaW = 2, betaW = 2)

# Hyperparameters for the inverse gamma (GAMMA^-1) prior of the hypervariances (tau^2)
# Default values determine the shape and scale of the prior distribution.
tau2 <- c(a1 = 5, a2 = 25)

# Set v_0, the ratio between the spike and slab variances (controls the spike strength)
# Default value sets a small spike.
gamma <- c(v0 = 0.00025)

# Hyperparameters for the inverse gamma (GAMMA^-1) prior for error variance (sigma^2)
# Specify the prior for error variance (relevant for Gaussian response).
# Default values control the prior distribution for error variance.
sigma2 <- c(b1 = 1e-4, b2 = 1e-4)

# Variance for the prior of xi
# This hyperparameter influences the prior for model parameters.
# No default value is specified in the provided information.

# Perform spike-and-slab regression on the data
sp_sl_model <- spikeAndSlab(y, X,family = "gaussian",
                            mcmc = list(nChains = 1, chainlength = 25000, burnin = 5000, thin = 1, sampleY = TRUE),
                            hyperparameters = list(w = w, tau2 = tau2, gamma = gamma, sigma2 = sigma2))

# Get the posterior means of the model parameters
coef <- sp_sl_model$postMeans

setwd("~/personal/spike-slab-analysis/reusable_functions/simulation_studies")
# Save the fitted model to a file for later use
save(sp_sl_model, file = file_name)



# Source the custom plot_function.R script
source("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis/reusable_functions/plot_function.R")


# # MCMC Chain Diagnostics
# sample <- sample(1:198, 10)  # Select 10 random samples from the MCMC chains
# 
# # Loop over the selected samples for diagnostic analysis
# for (s in sample) {
#   # Code for trace plots (currently commented out)
#   filename <- paste0("trace_plot_beta_", s, "(prior_beta(2,2)).png")
#   CairoPNG(filename = filename, width = 1800, height = 800,
#            pointsize = 12, bg = "white",  res = 200)
#   plot(sp_sl_model_R_beta_2_2$samples$beta[, s], main = bquote(beta[.(s)]))
#   dev.off()
#   # Assuming 'samples' is your MCMC sample matrix or data frame
#   mcmc_chain_list <- sp_sl_model_R_beta_2_2$samples$beta[, s]
#   
#   # Calculate the effective sample size for the MCMC chain
#   size = effectiveSize(mcmc_chain_list)
#   print(size)
#   
#   filename <- paste0("autocorrelation_plot_beta_", s, "prior_beta(2,2).png")
#   CairoPNG(filename = filename, width = 1800, height = 800,
#            pointsize = 12, bg = "white",  res = 200)
#   # Code for plotting the autocorrelation plot
#   autocorr.plot(mcmc_chain_list, main = bquote(beta[.(s)]))
#   
#   
#   # Add horizontal lines and labels to the autocorrelation plot
#   threshold = 0.1
#   abline(h = 0, col = "black", lwd = 2)
#   abline(h = threshold, col = "red", lty = 2)
#   abline(h = -threshold, col = "red", lty = 2)
#   text(25, 0.15, "0.1")
#   text(25, -0.15, "-0.1")
#   dev.off()
#   # Gelman diagnostic (commented out as it requires more than one chain)
#   # geweke_result <- geweke.diag(mcmc_chain_list)
#   # geweke.plot(mcmc_chain_list)
#   # print(geweke_result)
# }
# dev.off()
# graphics.off()
# # Note: Some code sections related to generating trace plots and Gelman diagnostics
# # are currently commented out or not provided. You may need to uncomment and adapt
# # those sections as needed for your analysis.
# 
# 
