X<-training_set[,-c(1:10)]
pca_result <- prcomp(X , scale = TRUE)
# variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
# cumulative_variance_explained <- cumsum(variance_explained)
pcs<-round(pca_result$rotation,2)

sig_features<-NULL

for(i in 1:50){
  pc<-pcs[,i]
  features<-names(pc[abs(pc)>0])
  columns = c(features,'R')
  temp_training_data<-training_set[,-c(2:10)]
  temp_training_data<-temp_training_data[columns]
  glm_model<-glm(R~.,data=temp_training_data)
  coef_summary<-as.data.frame(summary(glm_model)$coef)
  sig_features<-c(sig_features,rownames(coef_summary[coef_summary["Pr(>|t|)"]<0.05,]))
}
sig_features<-unique(sig_features)

temp_training_data<-training_set[,-c(1:10)]
temp_training_data<-temp_training_data[sig_features]
X_tmp<-temp_training_data[,-c(1:10)]
#choosing the phenotype
pheno='R'
y_tmp<-training_set[[pheno]]
#spike and slab regression on  the data
sp_sl_model_R_pca<-spikeAndSlab(y_tmp,X_tmp,mcmc=list(nChains=1,chainlength=25000,burnin=5000,thin=1,
                                              sampleY=TRUE))

#getting the posterior means of the parameters
coef<-sp_sl_model_R_pca$postMeans

temp_testing_data<-testing_set[,-c(1:10)]
temp_testing_data<-temp_testing_data[sig_features]
X_new_tmp<-as.matrix(temp_testing_data)
# pred<-X_new%*%beta
y_new_tmp<-testing_set[[pheno]]

#R_p calculation

model <- sp_sl_model_R_pca
pheno = "R"
beta <- model$postMeans$beta
# intercept<-beta[1]
# beta <- beta[-1]
pi <- c(1,model$postMeans$pV1)
coef <- data.frame(index=1:length(beta),beta,pi)
coef <- coef[order(coef$pi), ]


ss<- NULL
# c<-0.520236
for(prob in coef$pi){
  coef_new<-coef[coef$pi>prob,]
  
  df1_test <- X_new_tmp[ , coef_new$index]
  p_y_hat<-as.matrix(df1_test)%*%coef_new$beta
  ss_y_hat = 1-(sum((y_new_tmp-p_y_hat)^2)/n)/var(y_new_tmp)
  print(ss_y_hat)
  ss <- c(ss,ss_y_hat)
}

R_p <- ss
prob<-coef$pi
df_tmp<-data.frame(R_p,prob)
ridge<-ss_y_ridge[[pheno]]
lasso<-ss_y_lasso[[pheno]]
min_lim = min(c(ridge[,2],lasso[,2],ss))
max_lim = max(c(ridge[,2],lasso[,2],ss))
y_lim = c(min_lim,max_lim)
# y_lim = c(min(ss),max(max(ss)))
curve(custom_function(x, prob, ss), from = log(10^-4), to = log(max(prob)), n = 1000,
      xlab = expression(paste("log(Inclusion Probability) (", log(pi[~p]), ")")),
      ylab = expression(R[~p]^2), cex.lab = 0.7, cex.axis = 0.7,
      main = pheno, lwd = 3,ylim=y_lim,cex.main=1)
