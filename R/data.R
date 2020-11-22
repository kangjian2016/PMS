#'Data generation function
#'
#'A function used to generate data.
#'@export
#'@importFrom stats binomial
#'@importFrom stats poisson
#'@importFrom stats runif
#'@importFrom stats pnorm
#'@importFrom stats rnorm
#'@importFrom stats var
#'
#'@param n an integer represents the sample size.
#'@param p an integer which is the dimension of covariates.
#'@param rs signal to noise ratio.
#'@param beta_nonzero a vector corresponds to the values of non-zero coefficients in linear regression model.
#'@return a list of three variables, where "X" is a n by p matrix, n observations of p-dimensional covariates, "Y" is n observtions
#'        of the response.
#'@author Jie He \email{<jiehe@@umich.edu>}, Jian Kang \email{<jiankang@@umich.edu>}

data_gen = function(n,p,rs,beta_nonzero){
  x = matrix(rnorm(n*p),nrow = n, ncol = p)+matrix(rnorm(n, mean = 0,sd = sqrt(0.5)),nrow = n,ncol = p)
  Truebeta = rep(0,p)
  t_0 = length(beta_nonzero)
  Truebeta[1:t_0]=beta_nonzero
  Trueset=c(1:t_0)
  x=scale(x)
  y=x%*%Truebeta
  y=y+rnorm(n,0,sd=sqrt(var(y)/rs))
  return(list(x = x, y = y))
}



