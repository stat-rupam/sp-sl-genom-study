#code to generate simulated data
set.seed(123)
K = 1000 #number of SNPs
n = 5000 #sample size
theta = 0.2 #fraction  of causal SNPs
p = runif(1,0.05,0.5) #simulate random for binomial
beta = c(rnorm(K*theta),rep(0,K*(1-theta))) #generating betas
snps = matrix(rbinom(K * n, 2, p), nrow = 2*n, ncol = K) #generating the SNPs
error = rnorm(2*n) #generating the error term
Y = snps%*%beta+error #generating sample values

#calculate the actual R^2_p
Y_test = Y[(n+1):(2*n),]#separting the test set samples
R_p_act = 1 - (var(error[(n+1):(2*n)])/var(Y_test))

sim_data <- list()
sim_data[['Y']] <- Y
sim_data[['SNPs']] <- snps
sim_data[['R_p_act']] <- R_p_act
save(sim_data,file="sim_data_K_1000_theta_20.RData")