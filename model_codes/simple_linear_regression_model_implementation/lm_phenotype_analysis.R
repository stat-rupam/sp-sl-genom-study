library(dbplyr)
joined_table <- merge(phenod, cov, by = "IID",all.y = T)
joined_table <- merge(joined_table, geno, by = "IID",all.y = T)

col_means <- colMeans(joined_table[, -1], na.rm = TRUE)
df <- joined_table
df[, -1] <- lapply(df[, -1], function(x) ifelse(is.na(x), col_means[names(df[-1]) == names(x)], x))# Print the imputed data frame
genotypes<-colnames(geno)[-1]

phenotypes<-colnames(phenod)[-1]
covariates<-colnames(cov)[-1]

coef_df_list <- list()

df1=df[colnames(df)[-c(2:35)]]
df1['R'] = df['R']

for(pheno in phenotypes){
  coef_df =NULL
  for(gen in genotypes){
    df_name <- paste0(pheno,"_coef_df")
    formula=as.formula(paste0(pheno, " ~ ", paste(c(covariates,gen), collapse="+")))
    model <- lm(formula,df)
    # Print the model summary
    summary<-summary(model)
    print(pheno)
    print(summary)
    a=summary$coefficients[gen,]
    coef_df=rbind(coef_df,a)
  }
  coef_df=data.frame(coef_df,row.names = NULL)
  SNPs<-genotypes
  coef_df <- cbind(SNPs,coef_df)
  coef_df  <- coef_df[coef_df$Pr...t..<0.01,]
  coef_df_list[[df_name]] <- coef_df
}


#####
pheno="PC_RGB_3"
gen="rs28777"
genom2=c("rs5756492","rs1800407","rs1426654","rs12913832")
Interactions<-NULL
coef_df <- NULL
for (gen2 in genom2){
  df_name <- paste0(pheno,"_coef_df")
  formula=as.formula(paste0(pheno, " ~ ", paste(c(covariates, gen, gen2), collapse="+"), "+", gen, "*", gen2)
  )
  model <- lm(formula,df)
  # Print the model summary
  summary<-summary(model)
  Interactions<- c(Interactions,paste0(gen,":",gen2))
  a=summary$coefficients[paste0(gen,":",gen2),]
  coef_df=rbind(coef_df,a)
  print(a)
}
coef_df=data.frame(coef_df,row.names = NULL)
coef_df <- cbind(Interactions,coef_df)
coef_df  <- coef_df[coef_df$Pr...t..<0.01,]




