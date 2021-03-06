% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{data_gen}
\alias{data_gen}
\title{Data generation function}
\usage{
data_gen(n, p, rs, beta_nonzero)
}
\arguments{
\item{n}{an integer represents the sample size.}

\item{p}{an integer which is the dimension of covariates.}

\item{rs}{signal to noise ratio.}

\item{beta_nonzero}{a vector corresponds to the values of non-zero coefficients in linear regression model.}
}
\value{
a list of three variables, where "X" is a n by p matrix, n observations of p-dimensional covariates, "Y" is n observtions
       of the response.
}
\description{
A function used to generate data.
}
\author{
Jie He \email{<jiehe@umich.edu>}, Jian Kang \email{<jiankang@umich.edu>}
}
