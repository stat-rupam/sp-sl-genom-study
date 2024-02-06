library(BoomSpikeSlab)
library(coda)
library(dplyr)
library(spikeslab)
# Load training and testing data (temporary, for reusability; may be removed later)
# The paths to data files are specified here
load("~/personal/spike-slab-analysis/train_data_10000.RData")  # Load the training data
load("~/personal/spike-slab-analysis/test_data_10000.RData")   # Load the testing data

# Extract the predictors (features) and choose a specific phenotype for analysis
df1 <- training_set[,-c(1:10)]  # Extract the predictor variables
pheno = 'R'                  # Specify the phenotype of interest
df1[[pheno]] <- training_set[[pheno]]    # Extract the response variable
boom_sp_sl_R <- lm.spike(R ~. , df1,niter = 25000,
                         prior.information.weight = 10^-3,
                         diagonal.shrinkage = 0.00025)
coef_df <- summary(boom_sp_sl_R)$coef
coef_df <- data.frame(coef_df)
file = paste0("boom_spike_slab_output_",pheno,".RData")
save(boom_sp_sl_R, file=file)
              
#fit <- spikeslab(forumula,df1[-1], verbose=TRUE,n.iter1 = 500, n.iter2 = 100)
coef_df_list2 <- list()
for(pheno in phenotypes){
  df = coef_df_list[[pheno]]
  coef_df_list1[[pheno]] = df[1:10,]
}
# 
# #boom-spike-slab
# coef_list <- list()
# coef_list_lm <- list()
# 
# for (pheno in phenotypes){
#   formula=as.formula(paste0(pheno, " ~ ", "."))
#   df1=df[colnames(df)[-c(1:35)]]
#   df1[pheno] = df[pheno]
#   df1 <-as.data.frame(scale(df1))
#   model <- lm.spike(formula, data = df1,niter=100000,prior.information.weight = 10^-3,
#                     diagonal.shrinkage = 10^-5, prior.df = 3)
#   lm_model <- lm(formula,df1)
#   s= summary(model,burn=20000)
#   lm_s = summary(lm_model)
#   coef_list[[pheno]] = as.data.frame(s$coefficients)
#   coef_list_lm[[pheno]] = as.data.frame(lm_s$coefficients)
# }
# #save(coef_df_list1, file="boom_spike_slab_output.RData")
# 
# #coef_df_list2 <- list()
# 
# for(pheno in phenotypes){
#   df = coef_list[[pheno]]
#   names <- rownames(df)
#   rownames(df) <- NULL
#   data <- cbind(names,df)
#   data <-data[order(data$names), ]
#   data <- data[data$names %in% colnames(geno), ]
#   
#   df = coef_list_lm[[pheno]]
#   names <- rownames(df)
#   rownames(df) <- NULL
#   data1 <- cbind(names,df)
#   data1 <-data1[order(data1$names), ]
#   data1 <- data1[data1$names %in% colnames(geno), ]
# }
# 
# 
# #y=as.matrix(df1['R'],nrow=1)
# #X=as.matrix(df1[c(-1,-38)])
# #model <- spikeslab(y~X)
# 
# # Summary of the fitted model
# summary(model)
# # Extract MCMC samples
# mcmc_samples <- model$beta[,10]
# mcmc_object <- mcmc(as.matrix(mcmc_samples))
# traceplot(mcmc_object)
# # Access the MCMC samples
# print(mcmc_samples)
