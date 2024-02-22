library(rjags)
library(MCMCpack)

spike_slab_linear<-function(X, y, model_file, params, n_chain,chain_length,
                            burnin, thin=1, b1 = 1.0, b2 = 1.0, 
                            a_pi = 1.0, b_pi = 1.0, 
                            a1 = 5.0, a2 = 25.0, 
                            nu = 0.00025){
  # data creation
  N = dim(X)[1]
  p = dim(X)[2]
  
  dataList <- list(X = X, y = y, N = N, p = p, 
                   b1 = b1, b2 = b2, 
                   a_pi = a_pi, b_pi = b_pi, 
                   a1 = a1, a2 = a2, 
                   nu = nu)
  # Assuming y might be a dataframe column or in an incorrect shape
  dataList$y <- as.vector(dataList$y)  # Convert y to a vector if it's not already
  
  # Create a JAGS model
  model <- jags.model(textConnection(model_file), data = dataList, 
                      n.chains = n_chain)
  parallel.seeds("base::BaseRNG", n_chain)
  
  # Update the model (burn-in)
  update(model, burnin)
  
  # Sample from the posterior
  samples <- coda.samples(model, variable.names = params, 
                          n.iter = chain_length, thin = thin)
  retrun(samples)
}