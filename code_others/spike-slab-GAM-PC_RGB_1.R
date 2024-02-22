#install necessary libraries
install.packages(c("spikeSlabGAM", "coda", "Cairo", "parallel"));
# Load necessary libraries
library(spikeSlabGAM)  # Load the spikeSlabGAM package for the analysis
library(coda)          # Load the coda package for MCMC chain diagnostics
library(Cairo)         # Load the Cairo package for generating graphics
# Load the parallel package
library(parallel)

# Check the number of available CPU cores
num_cores <- detectCores()
print(num_cores)


# Load training and testing data (temporary, for reusability; may be removed later)
# The paths to data files are specified here
load("train_data_10000.RData")  # Load the training data
load("test_data_10000.RData")   # Load the testing data

# Extract the predictors (features) and choose a specific phenotype for analysis
X <- training_set[,-c(1:10)]  # Extract the predictor variables
pheno = 'PC_RGB_1'                  # Specify the phenotype of interest
y <- training_set[[pheno]]    # Extract the response variable

# Define hyperparameters for the spike-and-slab regression model
# Hyperparameters for the Beta-Beta prior for weights (w)
# Default values are for a uniform distribution.
w <- c(alphaW = 1, betaW = 1)

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
sp_sl_model_PC_RGB_1 <- spikeAndSlab(y, X, 
                              mcmc = list(nChains = 1, chainlength = 25000, burnin = 5000, thin = 1, sampleY = TRUE),
                              hyperparameters = list(w = w, tau2 = tau2, gamma = gamma, sigma2 = sigma2))

# Get the posterior means of the model parameters
coef <- sp_sl_model_PC_RGB_1$postMeans

file = paste0("sp_sl_model_",pheno,"_10000.RData")
# Save the fitted model to a file for later use
save(sp_sl_model_PC_RGB_1, file = file)

# MCMC Chain Diagnostics
sample <- sample(1:90, 10)  # Select 10 random samples from the MCMC chains

# Loop over the selected samples for diagnostic analysis
for (s in sample) {
  # Code for trace plots (currently commented out)
  # plot(sp_sl_model_PC_RGB_1$samples$beta[, s], main = bquote(beta[.(s)])
  
  # Assuming 'samples' is your MCMC sample matrix or data frame
  mcmc_chain_list <- sp_sl_model_PC_RGB_1$samples$beta[, s]
  
  # Calculate the effective sample size for the MCMC chain
  size = effectiveSize(mcmc_chain_list)
  print(size)
  
  # Code for plotting the autocorrelation plot
  autocorr.plot(mcmc_chain_list, main = bquote(beta[.(s)]))
  
  # Add horizontal lines and labels to the autocorrelation plot
  threshold = 0.1
  abline(h = 0, col = "black", lwd = 2)
  abline(h = threshold, col = "red", lty = 2)
  abline(h = -threshold, col = "red", lty = 2)
  text(25, 0.15, "0.1")
  text(25, -0.15, "-0.1")
  
  # Gelman diagnostic (commented out as it requires more than one chain)
  # geweke_result <- geweke.diag(mcmc_chain_list)
  # geweke.plot(mcmc_chain_list)
  # print(geweke_result)
}

# Note: Some code sections related to generating trace plots and Gelman diagnostics
# are currently commented out or not provided. You may need to uncomment and adapt
# those sections as needed for your analysis.




### Model File loading---not rquired to run
# for(pheno in phenotypes){
#   print(pheno)
#   file = paste("~/personal/spike-slab-analysis/sp_sl_model_",pheno,".RData")
#   file = gsub(" ", "", file)
#   load(file)
# }

