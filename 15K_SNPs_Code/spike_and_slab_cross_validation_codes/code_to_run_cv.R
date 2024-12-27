#code for model training and model performance check
# Load the sys library (assuming it's not a typo, sys is not a standard R package)
install.packages("sys","coda")
library(sys)
# Load necessary libraries
library(coda)          # Load the coda package for MCMC chain diagnostics
# library(Cairo)         # Load the Cairo package for generating graphics

#load required functions
source("~/personal/spike-slab-analysis/Cluster_codes/cross_validation_codes/cv_spike_slab_candella.R")
pheno <- "R"
n_fold <- 10
cv_opt_model <- cv_analysis_single_phenotype(data = training_set, 
                             phenotype = pheno, 
                             n_fold = n_fold, 
                             niter = niter, 
                             burn = burn, 
                             prior_information_weight = diagonal_shrinkage, 
                             diagonal_shrinkage = prior_information_weight)







# Source the custom plot_function.R script
#source("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis/reusable_functions/plot_function.R")

