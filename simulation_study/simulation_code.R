#code to generate simulated data
set.seed(123)
K = 50 #number of SNPs
n = 5000 #sample size
theta = 0.2 #fraction  of causal SNPs
p = runif(1,0.05,0.5) #simulate random for binomial
beta = c(rnorm(K*theta,sd=6),rep(0,K*(1-theta))) #generating betas
snps = matrix(rbinom(K * n, 2, p), nrow = 2*n, ncol = K) #generating the SNPs
error_sd = 11
error = rnorm(2*n,sd=error_sd) #generating the error term
Y = snps%*%beta+error #generating sample values

#calculate the actual R^2_p
Y_test = Y[(n+1):(2*n),] #separting the test set samples
R_p_act = 1 - (var(error)/var(Y))
print(var(Y))
sim_data <- list()
sim_data[['Y']] <- Y
sim_data[['SNPs']] <- snps
sim_data[['R_p_act']] <- R_p_act
sim_data[['p']] <- p
file_name_data <- paste0("sim_data_K_", K, "_theta_", (100 * theta), ".RData")
save(sim_data,file=file_name_data)

