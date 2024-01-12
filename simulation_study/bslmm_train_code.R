#code for model training and model performence check
# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)
# Load necessary libraries
library(monomvn)  # Load the spikeSlabGAM package for the analysis
library(coda)          # Load the coda package for MCMC chain diagnostics
library(Cairo)         # Load the Cairo package for generating graphics
library(doParallel)
library(hibayes)
cl <- makePSOCKcluster(50)
registerDoParallel(cl)
# Load training and testing data (temporary, for reusability; may be removed later)
# The paths to data files are specified here
K = 100
theta = 0.2
file = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
file_name = paste0("sim_data_K_",K,"_theta_",(100*theta),".RData")
folder = "simulated_data_v/"
file = paste0(file,folder,file_name)
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
Id<-1:length(y)
df<-data.frame(Id,y,X)
bslmm_model <- ibrm(y~.,data=df,M=X,M.id=1:length(y),method="BSLMM"
                      ,niter = 25000
                      ,nburn = 5000)
beta<-bslmm_model$beta
setwd("~/personal/spike-slab-analysis/reusable_functions/simulation_studies")
# Save the fitted model to a file for later use
file_name = paste0("bslmm_model_K_",K,"_theta_",(100*theta),".RData")
save(bslmm_model, file = file_name)
stopCluster(cl)