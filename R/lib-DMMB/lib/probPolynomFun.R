##===================================================================
# Simple function for nth order polynomial transformed through the  #
# probit function, so it is 'squashed' to be in the (0,1) interval  #
##===================================================================
probPolynomFun <- function(theta, X){
  deg <- NROW(theta)
  g   <- rep(0, NROW(X))
  for (i in 1:deg){
    g <- g + theta[i]*X^(deg-i)
  }
  Phi <- pnorm(g)
  return(Phi)
}