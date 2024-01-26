# Load the sys library (assuming it's not a typo, sys is not a standard R package)
library(sys)

# Source the custom plot_function.R script
source("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis/reusable_functions/plot_function.R")

# Source the ridge-lasso-sample-code.R script
source("C:/Users/rupam.basu/Documents/personal/spike-slab-analysis/reusable_functions/ridge-lasso-sample-code.R")

# Create a model object for further analysis
model <- sp_sl_model_R

# Specify the phenotype of interest
pheno = "R"

# Extract the estimated beta coefficients and pi values from the model
beta <- model$postMeans$beta
pi <- c(1, model$postMeans$pV1)
coef <- data.frame(index = 1:length(beta), beta, pi)
coef <- coef[order(coef$pi), ]

# Prepare the predictor variables for the testing set
X_new <- as.matrix(testing_set[,-c(1:10)])
y_new <- testing_set[[pheno]]

# Initialize an empty vector for storing results
ss <- NULL

# Loop through different probability thresholds
for (prob in coef$pi) {
  # Filter coefficients based on the current probability threshold
  coef_new <- coef[coef$pi > prob,]
  
  # Select the corresponding predictors
  df1_test <- X_new[, coef_new$index]
  
  # Calculate predicted values using selected coefficients
  p_y_hat <- as.matrix(df1_test) %*% coef_new$beta
  
  # Calculate R-squared for the current model
  ss_y_hat <- 1 - (sum((y_new - p_y_hat)^2) / n) / var(y_new)
  print(ss_y_hat)
  
  # Append the R-squared value to the results vector
  ss <- c(ss, ss_y_hat)
}

# Create a data frame to store R-squared values and associated probabilities
R_p <- ss
prob <- coef$pi
df <- data.frame(R_p, prob)

# Store the results in the 'sp_sl_results' list under the specified phenotype
sp_sl_results[[pheno]] <- df

# Retrieve Ridge and Lasso results (Befor running this do run ridge-lasso-sample-code.R)
ridge <- ss_y_ridge[[pheno]]
lasso <- ss_y_lasso[[pheno]]

# Determine the y-axis limits for the plot
min_lim = min(c(ridge[,2], lasso[,2], ss))
max_lim = max(c(ridge[,2], lasso[,2], ss))
y_lim = c(min_lim, max_lim)

# Define the output filename for the plot
filename <- paste0("sp_sl_R2_plots_", pheno, ".png")

# Generate a PNG plot using Cairo
CairoPNG(filename = filename, width = 1000, height = 1000,
         pointsize = 12, bg = "white", res = 200)

# Create a custom plot using the R_p_plot function
curve(R_p_plot(x, prob, ss), from = log(10^-4), to = log(max(prob)), n = 1000,
      xlab = expression(paste("log(Inclusion Probability) (", log(pi[~p]), ")")),
      ylab = expression(R[~p]^2), cex.lab = 0.7, cex.axis = 0.7,
      main = pheno, lwd = 3, ylim = y_lim, cex.main = 1)

# Optionally, save the 'sp_sl_results' to a file or close the graphics device
# save(sp_sl_results, file = "sp_sl_model_performance.RData")
# graphics.off()
