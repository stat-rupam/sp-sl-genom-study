library(BoomSpikeSlab)
library(coda)
library(dplyr)
library(spikeslab)
joined_table <- merge(phenod, cov, by = "IID",all.y = T)
joined_table <- merge(joined_table, geno, by = "IID",all.y = T)

col_means <- colMeans(joined_table[, -1], na.rm = TRUE)
df <- joined_table

# Fill missing values with column means
df[colSums(is.na(df)) > 0] <- lapply(df[colSums(is.na(df)) > 0], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})


#df[, -1] <- lapply(df[, -1], function(x) ifelse(is.na(x), col_means[names(df[-1]) == names(x)], x))# Print the imputed data frame
phenotypes<-colnames(phenod)[-1]
coef_df_list <- list()

pheno = "limbalring1"
df1=df[colnames(df)[-c(1:35)]]
df1[pheno] = df[pheno]
fit <- spikeslab(a_lab ~. , df1[-1], verbose=TRUE,n.iter1 = 5000, n.iter2 = 1000)
coef_df_list[[pheno]] <- fit$summary[1:2]

save(coef_df_list, file="spike_slab_output.RData")

#fit <- spikeslab(forumula,df1[-1], verbose=TRUE,n.iter1 = 500, n.iter2 = 100)
coef_df_list2 <- list()
for(pheno in phenotypes){
  df = coef_df_list[[pheno]]
  coef_df_list1[[pheno]] = df[1:10,]
}

#boom-spike-slab
coef_list <- list()
coef_list_lm <- list()

for (pheno in phenotypes){
  formula=as.formula(paste0(pheno, " ~ ", "."))
  df1=df[colnames(df)[-c(1:35)]]
  df1[pheno] = df[pheno]
  df1 <-as.data.frame(scale(df1))
  model <- lm.spike(formula, data = df1,niter=100000,prior.information.weight = 10^-3,
                    diagonal.shrinkage = 10^-5, prior.df = 3)
  lm_model <- lm(formula,df1)
  s= summary(model,burn=20000)
  lm_s = summary(lm_model)
  coef_list[[pheno]] = as.data.frame(s$coefficients)
  coef_list_lm[[pheno]] = as.data.frame(lm_s$coefficients)
}
#save(coef_df_list1, file="boom_spike_slab_output.RData")

#coef_df_list2 <- list()

for(pheno in phenotypes){
  df = coef_list[[pheno]]
  names <- rownames(df)
  rownames(df) <- NULL
  data <- cbind(names,df)
  data <-data[order(data$names), ]
  data <- data[data$names %in% colnames(geno), ]
  
  df = coef_list_lm[[pheno]]
  names <- rownames(df)
  rownames(df) <- NULL
  data1 <- cbind(names,df)
  data1 <-data1[order(data1$names), ]
  data1 <- data1[data1$names %in% colnames(geno), ]
}


#y=as.matrix(df1['R'],nrow=1)
#X=as.matrix(df1[c(-1,-38)])
#model <- spikeslab(y~X)

# Summary of the fitted model
summary(model)
# Extract MCMC samples
mcmc_samples <- model$beta[,10]
mcmc_object <- mcmc(as.matrix(mcmc_samples))
traceplot(mcmc_object)
# Access the MCMC samples
print(mcmc_samples)
