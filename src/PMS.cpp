#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//'@importFrom Rcpp evalCpp
//'@useDynLib PMS, .registration=TRUE

//'@export
// [[Rcpp::export]]
Rcpp::List fast_PMS_cpp(arma::mat &x, arma::vec& y, double theta=1e-3){
  arma::mat XXt = x*x.t();
  XXt.diag() += theta;
  arma::vec screen_stat = arma::solve(XXt,y,arma::solve_opts::likely_sympd);
  screen_stat = x.t()*screen_stat;
  arma::vec abs_screen_stat = arma::abs(screen_stat);
  arma::uvec od = arma::sort_index(abs_screen_stat,"descend");
  return Rcpp::List::create(Named("pms_select") = od+1L,
                            Named("pms_stat") = screen_stat);
}

//'@export
// [[Rcpp::export]]
Rcpp::List fast_data_gen(int n, int p, double R2, arma::vec& beta_nonzero){
  arma::mat x = arma::randn<arma::mat>(n,p);
  arma::vec z = arma::randn<arma::vec>(n);
  for(int i=0;i<x.n_cols;i++){
  	x.col(i) += z;
  }
  arma::vec y = x.cols(0,beta_nonzero.n_elem-1L)*beta_nonzero;
  double var_y = arma::var(y);
  y += arma::randn(n)*sqrt(var_y*(1.0 - R2)/R2);

  return Rcpp::List::create(Named("x") = x,
                            Named("y") = y);
}
