library(rjags)
library(MCMCpack)
source("C:\\Users\\Rupam Basu\\Desktop\\spike_and_slab_jags_code.R")
modelString <- "
model {
  # Priors
  
  # sigma^2 prior
  sigma2 ~ dgamma(b1, b2)
  
  # prior for pi
  pi ~ dbeta(a_pi, b_pi)
  
  # prior for beta
  for (i in 1:p) {
    # prior for delta
    delta[i] ~ dbern(pi)
    
    #prior for tau
    tau[i] ~ dgamma(a1, a2)
    # Conditional distribution of beta[i] depending on delta[i]
    beta[i] ~ dnorm(0.0, pow(nu * (1/tau[i]) * (1 - delta[i]) + delta[i] * (1/tau[i]), 0.5))
  }
  
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], pow(sigma2, -0.5))
    mu[i] <- inprod(X[i,], beta)
  }
}
"

# Parameters to monitor
params <- c("beta", "sigma2", "pi", "delta", "tau")

model<-spike_slab_linear(X,y,modelString,params,
                         3,5000,1000)
# # Create a JAGS model
# model <- jags.model(textConnection(modelString), data = dataList, 
#                     n.chains = 3)
# 
# # Update the model (burn-in)
# update(model, 1000)
# 
# # Sample from the posterior
# samples <- coda.samples(model, variable.names = params, n.iter = 5000)

# Check results
s<-summary(samples)$statistics
View(data.frame(true_beta,est_beta=s[grep("beta",rownames(s)),1],
                inc_prob=1-s[grep("delta",rownames(s)),1]))
plot(samples)
