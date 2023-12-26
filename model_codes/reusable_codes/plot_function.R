# Load the Cairo library
library(Cairo)

# Explanation for the 'Cairo' library:
# Cairo is a graphics device for R that provides high-quality graphics output, including support for various file formats and antialiasing.

# Define the 'R_p_plot' function
R_p_plot <- function(x, prob, response) {
  output = NULL
  n <- length(prob)  # Get the length of the probability vector
  
  for (j in x) {  # Loop through the values of 'x'
    for (i in 1:n) {  # Loop through the indices of 'prob'
      if (j <= log(prob[i])) {  # Check if the current 'x' is less than or equal to the log probability
        output = c(output, response[i])  # Append the corresponding response value to 'output'
        break  # Exit the inner loop once a match is found
      }
    }
  }
  return(output)  # Return the 'output' vector containing response values
}
# Define the 'R_p_plot' function
R_p_plot_by_size <- function(x, size, response) {
  output = NULL
  n <- length(size)  # Get the length of the probability vector
  
  for (j in x) {  # Loop through the values of 'x'
    for (i in 1:n) {  # Loop through the indices of 'prob'
      if (j < size) {  # Check if the current 'x' is less than or equal to the lsize
        output = c(output, response[i])  # Append the corresponding response value to 'output'
        break  # Exit the inner loop once a match is found
      }
    }
  }
  return(output)  # Return the 'output' vector containing response values
}
# Explanation for the 'R_p_plot' function:
# This function takes three arguments: 'x', 'prob', and 'response'.
# - 'x' is a vector of values.
# - 'prob' is a vector of probabilities.
# - 'response' is a vector of response values associated with the probabilities.
# The function iterates through 'x' and for each value, it looks up the corresponding probability and appends the associated response value to the 'output' vector.
# The function returns the 'output' vector, which can be used for plotting.
