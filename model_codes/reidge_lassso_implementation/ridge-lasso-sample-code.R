library(glmnet)
library(Cairo)
#declaring two list for saving the results for each phenotypes
ss_y_ridge<-list()
ss_y_lasso<-list()
# The paths to data files are specified here
load("~/personal/spike-slab-analysis/train_data_2000.RData")  # Load the training data
load("~/personal/spike-slab-analysis/test_data_2000.RData")   # Load the testing data
#name of the phenotypes
phenotypes<-c("R","G","B","C","PC_RGB_1","PC_RGB_2","PC_RGB_3","L_lab","a_lab","b_lab")

#ridge regression and lasso regression fitting
for(pheno in phenotypes){
  df1 <- training_set  # Use the entire training 
  # Scale only the predictor variables, not the response variable
  df1 <- scale(df1[, -c(1:34)])
  # pheno = 'R'
  #df1[pheno] <- training_set[pheno]
  
  X=as.matrix(df1)
  y = as.matrix(training_set[pheno])
  
  # Set the alpha parameter to 0 for Ridge regression (L2 penalty)
  lambda_values <- seq(0.01, 100, length = 100)
  ridge_model <- glmnet(X, y, alpha = 0, family = "gaussian")
  lasso_model <- glmnet(X, y, alpha = 1, family = "gaussian")
  
  df2 <- testing_set  # Use the entire training set
  
  # Scale only the predictor variables, not the response variable
  df2 <- scale(df2[, -c(1:34)])
  #df1[pheno] <- training_set[pheno]
  X_new = as.matrix(df2)
  y_new <- testing_set[pheno]
  n <-dim(y)[1]
  #extracting the lambdas
  lambda_ridge<- ridge_model$lambda
  lambda_lasso<- lasso_model$lambda
  
  #prediction for different choices of lambdas
  predictions_ridge <- predict(ridge_model, newx = X_new, type="response")
  predictions_lasso <- predict(lasso_model, newx = X_new, type="response")
  
  #calculating R^2 for different lambdas for ridge
  R2_ridge<-as.vector(apply(predictions_ridge, 2, function(col) 1-sum((col - y_new)^2)/n/var(y_new)))
  ss_y_ridge[[pheno]]<-data.frame(lambda_ridge,R2_ridge)
  
  #calculating R^2 for different lambdas for lasso
  R2_lasso<-as.vector(apply(predictions_lasso, 2, function(col) 1-sum((col - y_new)^2)/n/var(y_new)))
  ss_y_lasso[[pheno]]<-data.frame(lambda_lasso,R2_lasso)
}

#This part is for plotting the the R^2 vs lambda values
par(mfrow=c(1,2),mar = c(5, 5,0.8,0.8))
for(pheno in phenotypes){
  # mar = c(1.5, 1.7,0.8,0.8)
  print(pheno)
  ridge<-ss_y_ridge[[pheno]]
  lasso<-ss_y_lasso[[pheno]]
  # ss_ridge<-c(ss_ridge,max(ridge[,2]))
  # ss_lasso<-c(ss_lasso,max(lasso[,2]))
  # lambda_ridge<-c(log(lambda_ridge),ridge[ridge[,2]==max(ridge[,2]),1])
  # lambda_lasso<-c(log(lambda_lasso),lasso[lasso[,2]==max(lasso[,2]),1])
  min_lim = min(c(ridge[,2],lasso[,2]))
  max_lim = max(c(ridge[,2],lasso[,2]))
  y_lim = c(min_lim,max_lim)
  filename <- paste0("ridge_lasso_R2_plots_", pheno, ".png")
  CairoPNG(filename = filename, width = 1000, height = 1000,
           pointsize = 12, bg = "white",  res = 200)
  par(mfrow=c(1,2),mar = c(5, 5,0.8,0.8))
  plot(log(ridge[,1]),ridge[,2],type='l',main=paste(pheno,"-Ridge"),
       xlab = expression(log(lambda)),
       ylab = expression(R[~p]^2),lwd=3,col='blue',cex.main=0.8,ylim=y_lim)
  plot(log(lasso[,1]),lasso[,2],type='l',main=paste(pheno,"-Lasso"),
       xlab = expression(log(lambda)),
       ylab = expression(R[~p]^2),lwd=3,col='green',cex.main=0.8,ylim=y_lim)
}
graphics.off()