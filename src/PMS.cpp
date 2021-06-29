#include <RcppArmadillo.h>
#ifndef NDEBUG
#define NDEBUG
#endif
#include "RNifti.h"
#include "RNiftiAPI.h"



// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



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
  for(arma::uword i=0;i<x.n_cols;i++){
  	x.col(i) += z;
  }
  arma::vec y = x.cols(0,beta_nonzero.n_elem-1L)*beta_nonzero;
  double var_y = arma::var(y);
  y += arma::randn(n)*sqrt(var_y*(1.0 - R2)/R2);
  return Rcpp::List::create(Named("x") = x,
                            Named("y") = y);
}

//'@export
// [[Rcpp::export]]
Rcpp::List fast_PMS_local_spatial(arma::mat& x, arma::vec& y,
                                  arma::mat& coords,
                                  arma::mat& neighbors,
                                  arma::uvec& num_neighbors,
                                  arma::vec& rho, double theta=1e-3){

  //create Lambda*t(X)
  //exp(-rho(v - v')^2)

  arma::mat lx = arma::zeros<arma::mat>(x.n_rows,x.n_cols);

  arma::mat XXt = x*lx.t();
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
arma::mat fast_read_imgs_mask(CharacterVector& imgfiles,
                              CharacterVector& maskfile,
                              int verbose = 100)
{
  const RNifti::NiftiImage mask(maskfile);
  arma::uvec mask_uvec = Rcpp::as<arma::uvec>(mask.toArray());
  arma::uvec mask_idx = arma::find(mask_uvec>0);

  int num_imgs = imgfiles.length();
  int num_voxels = mask_idx.n_elem;
  if(verbose < 1){
  	verbose = 1;
  }

  arma::mat img_mat(num_voxels,num_imgs);
  for(int id=0;id<imgfiles.length();id++){
    const RNifti::NiftiImage image(Rcpp::as<std::string>(imgfiles[id]));
    arma::vec img_vec = Rcpp::as<arma::vec>(image.toArray());
    img_mat.col(id)= img_vec.elem(mask_idx);
    if((id+1)%verbose==0){
    	std::cout << id+1 << std::endl;
    }
  }

  return img_mat;
}




