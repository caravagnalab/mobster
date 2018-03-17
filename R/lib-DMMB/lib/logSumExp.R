logSumExp <- function(x) {
  ##=========================================================
  # Function to compute the log sum exp trick for avoiding  #
  # numeric underflow and have numeric stability in the     #
  # computations of really small numbers.                   #
  ##=========================================================
  # Computes log(sum(exp(x))
  offset <- max(x)
  return(log(sum(exp(x - offset))) + offset)
}