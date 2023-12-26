library(genio)
library(caret)
library(glmnet)
#change the working directory
setwd("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis")

#load the covariates and phenotype values
load("~/personal/spike-slab-analysis/pigmentation_data_for_Sabya.RData")

#load the genotype values
data_genio = read_plink("cand16pigm2_geno")
format(object.size(data_genio), units = "auto")
genomat=t(data_genio$X); # 6237 rows and 127857 cols
bim=data_genio$bim;
rm(data_genio)
k=200;
f=sample.int(nrow(bim), size=k, replace = FALSE)
subgeno=data.frame(genomat[,f]);
subgeno$IID=rownames(subgeno)

#join  the covariates, phenotypes and the genotypes tables
joined_table <- merge(phenod, cov, by = "IID",all.y = T)
joined_table <- merge(joined_table, subgeno, by = "IID",all.y = T)

col_means <- colMeans(joined_table[, -1], na.rm = TRUE)
df <- joined_table

# Fill missing values with column means
df[colSums(is.na(df)) > 0] <- lapply(df[colSums(is.na(df)) > 0], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})


#df[, -1] <- lapply(df[, -1], function(x) ifelse(is.na(x), col_means[names(df[-1]) == names(x)], x))# Print the imputed data frame
phenotypes<-colnames(phenod)[-1]
df_modified<-drop_highly_correlated(df[, -c(1:47)],0.2)
# df_modified["IID"]<-df["IID"]

#name of the phenotypes
phenotypes<-c("R","G","B","C","PC_RGB_1","PC_RGB_2","PC_RGB_3","L_lab","a_lab","b_lab")
df_modified<-scale(df_modified[-1])
df_modified<-as.data.frame(df_modified)
df_modified["IID"]<-df["IID"]

df_modified1<-merge(df[c("IID",phenotypes)], df_modified, by = "IID",all.y = T)
# df_modified1[,colnames(df[,c(36:47)])]<-df[,c(36:47)]
df_modified1<-df_modified1[,-1]
#train-test split
set.seed(123)  # Set a seed for reproducibility
total_rows <- nrow(df_modified1)

test_size <- 0.2
total_rows<-dim(df)[1]
num_test_rows <- round(test_size * total_rows)
test_indices <- sample(total_rows, num_test_rows)
training_set <- df_modified1[-test_indices, ]
testing_set <- df_modified1[test_indices, ]
# save(training_set,file="train_data_2000.RData")
# save(testing_set,file="test_data_2000.RData")
