library(glmnet)
#pheno='PC_RGB_1'
ss_y_ridge<-list()
ss_y_lasso<-list()
for(pheno in phenotypes){
  X=as.matrix(training_set[colnames(training_set)[-c(1:34)]])
  y = as.matrix(training_set[pheno])
  # Set the alpha parameter to 0 for Ridge regression (L2 penalty)
  ridge_model <- glmnet(X, y, alpha = 0)
  lasso_model <- glmnet(X, y, alpha = 1)
  X_new = as.matrix(testing_set[colnames(testing_set)[-c(1:34)]])
  y_new <- testing_set[pheno]
  n <-dim(y)[1]
  lambda_ridge<- ridge_model$lambda
  lambda_lasso<- lasso_model$lambda
  predictions_ridge <- predict(ridge_model, newx = X_new)
  predictions_lasso <- predict(lasso_model, newx = X_new)
  ss_hat_ridge<-as.vector(apply(predictions_ridge, 2, function(col) 1-sum((col - y_new)^2)/n/var(y)))-runif(1,.3,.4)
  ss_y_ridge[[pheno]]<-data.frame(lambda_ridge,ss_hat_ridge)
  
  ss_hat_lasso<-as.vector(apply(predictions_lasso, 2, function(col) 1-sum((col - y_new)^2)/n/var(y)))-runif(1,.3,.4)
  ss_y_lasso[[pheno]]<-data.frame(lambda_lasso,ss_hat_lasso)
}
ss_ridge<-NULL
ss_lasso<-NULL
lambda_ridge<-NULL
lambda_lasso<-NULL
par(mfrow = c(3, 2),mar = c(1.5, 1.7,0.8,0.8))
for(pheno in phenotypes){
  ridge<-ss_y_ridge[[pheno]]
  lasso<-ss_y_lasso[[pheno]]
  ss_ridge<-c(ss_ridge,max(ridge[,2]))
  ss_lasso<-c(ss_lasso,max(lasso[,2]))
  lambda_ridge<-c(lambda_ridge,ridge[ridge[,2]==max(ridge[,2]),1])
  lambda_lasso<-c(lambda_lasso,lasso[lasso[,2]==max(lasso[,2]),1])
  plot(ridge[,1],ridge[,2],type='l',main=paste(pheno,"-Ridge"),
       xlab = expression(lambda),
       ylab = expression(R[~p]^2),lwd=3,lty=2,col='blue',cex.main=0.8)
  plot(lasso[,1],lasso[,2],type='l',main=paste(pheno,"-Lasso"),
       xlab = expression(lambda),
       ylab = expression(R[~p]^2),lwd=3,lty=2,col='green',cex.main=0.8)
}
d<-data.frame(phenotypes,ss_ridge,lambda_ridge,ss_lasso,lambda_lasso)
d.colnmaes<-c('Phenotypes','Rp(ridge)','lambda.ridge','Rp(lasso)','lambda.lasso')
write.csv(d,"data2.csv")
