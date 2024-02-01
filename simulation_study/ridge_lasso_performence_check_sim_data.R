library(glmnet)
library(Cairo)
# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)

K = 1000
theta = 0.2
n = 5000
folder_data_index = "iii" #iii,iv,v
folder_model_index = 70 #70,50,30
non_zero_coef_index = 1:(K*theta)

file = "~/personal/spike-slab-analysis/reusable_functions/simulation_studies/"
folder_data = paste0("simulated_data_" , folder_data_index , "/")
file_name_data = paste0("sim_data_K_",K,"_theta_",(100*theta),".RData")
file_name_model = paste0("sp_sl_model_K_",K,"_theta_",(100*theta),".RData")
folder_model = paste0("model_output_R_p_" , folder_model_index , "/")
file_name_data = paste0(file,folder_data,file_name_data)
file_name_model = paste0(file,folder_model,file_name_model)
load(file_name_model)
load(file_name_data)

#actual R_p
R_p_act = sim_data[['R_p_act']]
SNPS = scale(sim_data[['SNPs']])
training_set = SNPS[1:n,]
testing_set = SNPS[(n+1):(2*n),]
y_train = sim_data[['Y']][1:n]
y_test = sim_data[['Y']][(n+1):(2*n)]

#ridge regression and lasso regression fitting
df1 <- training_set  # Use the entire training 
# Scale only the predictor variables, not the response variable
X = as.matrix(df1)
y = as.matrix(y_train)

# Set the alpha parameter to 0 for Ridge regression (L2 penalty)
lambda_values <- seq(0.01, 100, length = 100)
ridge_model <- glmnet(X, y, alpha = 0, family = "gaussian")
lasso_model <- glmnet(X, y, alpha = 1, family = "gaussian")

df2 <- testing_set  # Use the entire training set

X_new = as.matrix(df2)
y_new <- y_test

#extracting the lambdas
lambda_ridge<- ridge_model$lambda
lambda_lasso<- lasso_model$lambda

#prediction for different choices of lambdas
predictions_ridge <- predict(ridge_model, newx = X_new, type="response")
predictions_lasso <- predict(lasso_model, newx = X_new, type="response")

#calculating R^2 for different lambdas for ridge
R2_ridge<-as.vector(apply(predictions_ridge, 2, function(col) 1-sum((col - y_new)^2)/n/var(y_new)))
# ss_y_ridge[[pheno]]<-data.frame(lambda_ridge,R2_ridge)

#calculating R^2 for different lambdas for lasso
R2_lasso<-as.vector(apply(predictions_lasso, 2, function(col) 1-sum((col - y_new)^2)/n/var(y_new)))
# ss_y_lasso[[pheno]]<-data.frame(lambda_lasso,R2_lasso)
non_zero_common_lasso <- apply(coef(lasso_model), 2, function(x) length(intersect(non_zero_coef_index,
                                                                      which(x != 0))))
num_var_lasso = apply(coef(lasso_model), 2, function(x) sum(x != 0))

print(data.frame(R_p_act, max(R2_ridge), max(R2_lasso), 
                 number_of_var=num_var_lasso[R2_lasso==max(R2_lasso)],
                 common = non_zero_common_lasso[R2_lasso==max(R2_lasso)]
))
#print(R2_lasso[num_var_lasso<=K*theta])
