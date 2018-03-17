##=======================================================================
# bin2ndProbRegLik - Binomial distributed Probit Regression Likelihood  #
# function. The expression for the log likelihood is:                   #
#   log L(theta;x) = Sum_{l=1}^{L}log(Bin(n_{l}, Phi(g(x_{l});theta)))  #
#                                                                       #
# where g(x)  = ax^2 + bx + c                                           #
#       theta = (a,b,c)                                                 #
#       n     = (m, t) --> m is the number of successes in t trials     #
#                                                                       #
# Two modes are provided, for computing the log likelihoods and the     #
# derivatives wrt to the parameters of the polynomials                  #
# and for computing the approximate predictive probability.             #
#                                                                       #
# Usage: y          <- bin2ndProbRegLik(theta, D, mode=1)               #
#        (da,db,dc) <- bin2ndProbRegLik(theta, D, mode=2)               #
#                                                                       #
# Input:                                                                #
#     theta   is a vector of the parameters of the function             #
#     D       is a 3 by L dimensional matrix containing:                #
#         X   is a vector (of length L) with the corresponding values   #
#               at each position                                        #
#         t   is a vector (of length L) with the number of trials, in   #
#               the context of Binomal distributed data                 #
#         m   is a vector (of length L) with the number of successes in #
#               the t trials                                            #
# Output:                                                               #
#     y       is the log likelihood                                     #
#                                                                       #
#     da      is the derivative wrt to the a parameter                  #
#     db      is the derivative wrt to the b parameter                  #
#     dc      is the derivative wrt to the c parameter                  #
#                                                                       #
##=======================================================================
bin2ndProbRegLik <- function(theta, D, mode=1){
  X   <- D[1,]            # Data of length L
  t   <- D[2,]            # Number of trials for the corresponding X
  m   <- D[3,]            # Number of successes for the corresponding t
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  Phi <- pnorm(g)         # Squash the function g() to the (0,1) interval
  
  if (mode==1){           # Compute the log likelihood
    res <- sum(dbinom(x=m, size=t, prob=Phi, log=TRUE))
    return(res)
  }else if (mode==2){     # Compute derivatives wrt to the parameters
    N   <- dnorm(g)       # Density function of the Normal distribution
    da  <- sum(N*X^2 * (m-t*Phi)/(Phi*(1-Phi))) # Der wrt to a parameter
    db  <- sum(N*X * (m-t*Phi)/(Phi*(1-Phi)))   # Der wrt to b parameter
    dc  <- sum(N * (m-t*Phi)/(Phi*(1-Phi)))     # Der wrt to c parameter
    return(c(da, db, dc))
  }
}
