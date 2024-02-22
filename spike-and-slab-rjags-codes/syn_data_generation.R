set.seed(123) # For reproducibility

N <- 2000 # Number of observations
p <- 1000  # Number of predictors

# Generate predictors
X <- matrix(rnorm(N * p), ncol = p)

# True beta coefficients - with some set to 0 to simulate irrelevant predictors
true_beta <- sample(1:5,200,replace = T)
true_beta <- c(true_beta,rep(0,800))

# Generate noise
sigma <- 1
noise <- rnorm(N, mean = 0, sd = sigma)

# Generate outcome variable
y <- X %*% true_beta + noise

# Data list for RJAGS
dataList <- list(X = X, y = y, N = N, p = p, 
                 b1 = 1.0, b2 = 1.0, 
                 a_pi = 1.0, b_pi = 1.0, 
                 a1 = 5.0, a2 = 25.0, 
                 nu = 0.00025)
# Assuming y might be a dataframe column or in an incorrect shape
dataList$y <- as.vector(dataList$y)  # Convert y to a vector if it's not already

# Verify dimensions
dataList$N <- length(dataList$y)  # Ensure N matches the length of y
dataList$p <- ncol(dataList$X)  # Ensure p matches the number of predictors in X
