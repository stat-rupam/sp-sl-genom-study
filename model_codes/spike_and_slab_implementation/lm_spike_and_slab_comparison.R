library(caret)
set.seed(123)  # Set a seed for reproducibility
total_rows <- nrow(df)
df<-scale(df[-1])
test_size <- 0.2
num_test_rows <- round(test_size * total_rows)
test_indices <- sample(total_rows, num_test_rows)
training_set <- df[-test_indices, ]
testing_set <- df[test_indices, ]

ss_y<-list()
for(pheno in phenotypes){
  #pheno='R'
  formula=as.formula(paste0(pheno, " ~ ", "."))
  df1=training_set[colnames(training_set)[-c(1:34)]]
  df1[pheno] = training_set[pheno]
  #df1 <-as.data.frame(df1)
  model <- lm.spike(formula, data = df1,niter=100000,prior.information.weight = 10^-3,
                    diagonal.shrinkage = 10^-5, prior.df = 3)
  df1_test = testing_set[colnames(testing_set)[-c(1:34)]]
  y <- testing_set[pheno]
  n <-dim(y)[1]
  #df1_test <-as.data.frame(scale(df1_test))
  # pred<-predict(model,newdata=df1_test,burn=20000,means.only=TRUE)
  # y_hat <- pred[,80000]
  #n<-length(y)
  #ss_y_hat<- 1- sum((y-y_hat)^2)/n
  s= summary(model,burn=20000)
  
  coef<-as.data.frame(s$coefficients)
  names <- rownames(coef)
  rownames(coef) <- NULL
  coef <- cbind(names,coef)
  intercept <- coef[match('(Intercept)', coef$names), 4]
  coef <- coef[-match('(Intercept)', coef$names),]
  inc_probs = rev(unique(coef[,6]))
  coef_name<-coef$names
  ss <- NULL
  for(prob in inc_probs[c(-1,-length(inc_probs))]){
    coef_new<-coef[coef[,6]>prob,]
    # df1_test["(Intercept)"]<-rep(1,n)
    df1_test <- df1_test[ , c(names(df1_test)[names(df1_test) != "(Intercept)"])]
    
    X<-df1_test[,match(coef_new$names,colnames(df1_test))]
    #View(X)
    coef_new <-na.omit(coef[match(coef_name, coef_new$names), 4])
    
    p_y_hat<-intercept + as.matrix(X)%*%as.matrix(coef_new)
    ss_y_hat = 1-(sum((y-p_y_hat)^2)/n)/var(y)
    print(ss_y_hat)
    ss <- c(ss,ss_y_hat)
  }
  p_y_hat<-intercept + as.matrix(X)%*%as.matrix(coef_new)
  ss_y_hat = 1-(sum((y-p_y_hat)^2)/n)/var(y)
  ss <- c(ss,ss_y_hat)
  ss_y[[pheno]] <- data.frame(ss, inc_probs[-1])
  lm_model <- lm(formula,df1)
  coef_name<- rownames(as.data.frame(summary(lm_model)$coef))
  pred1 <- predict(lm_model,df1_test)
  ss_y_hat<- 1- (sum((y-pred1)^2)/n)/var(y)
  ss_y[[paste(pheno,'_lm')]] <- ss_y_hat
}



lm_model <- lm(PC_RGB_1~.,df1)
coef_name<- rownames(as.data.frame(summary(lm_model)$coef))
pred1 <- predict(lm_model,df1_test)
ss_y_hat<- 1- sum((y-pred1)^2)/n


#plots

custom_function <- function(x, prob, response) {
  output = NULL
  n <- length(prob)
  for (j in x){
    for (i in 1:n) {
      if (j <= prob[i]) {
        output = c(output,response[i])
        break
      }
    }
  }
  return(output)
}
data<-read.csv("data.csv")
i=1
#par(mfrow = c(2, 2),mar = c(1.5, 1.7,0.8,0.8))
for(pheno in phenotypes){
  name = paste(pheno,".jpeg")
  jpeg(name)
  d<-ss_y[[pheno]]
  ss_y_hat<-as.vector(d['ss'])$ss
  prob<-as.vector(d["inc_probs..1."])$inc_probs..1.
  y_lim = c(min(ss_y_hat),max(max(ss_y_hat),data[i,'R_p.lm.']))
  curve(custom_function(x, prob, ss_y_hat), from = 0, to = max(prob), n = 1000,
        xlab = expression(paste("Inclusion Probability (", pi[~p], ")")),
        ylab = expression(R[~p]^2), cex.lab = 0.7, cex.axis = 0.7,
        main = pheno, lwd = 3,ylim=y_lim,cex.main=1)
  
  abline(h = data[i,'R_p.lm.'], col = "red", lty = 2, lwd = 3)
  
  # Adding the legend
  legend("topright", legend = expression(paste(R[~p]^2, " For multiple linear regresssion model"))
         , col = "red", lty = 2, lwd = 3, cex=0.8)
  i=i+1
  
  dev.off()
}

