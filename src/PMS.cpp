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

//'Fast posterior mean screening without prior knowledge
//'@param x: n x p covariate matrix
//'@param y: n x 1 response vector
//'@param theta: positive number control the prior precision. Default value is 1e-3.
//'@return a list of two variables,
//'\describe{
//' \item{pms_select}{indices of the predictors from most important to least important}
//' \item{pms_stat}{a p by 1 vector of posterior mean statistics}
//'}
//'@author Jian Kang \email{<jiankang@@umich.edu>}
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


//'Simulate data from linear regression model
//'@param n an integer represents the sample size.
//'@param p an integer which is the dimension of covariates.
//'@param R2 a number between 0 and 1 to indicate the signal to noise ratio.
//'@param beta_nonzero a vector of q (q > 0) non-zero coefficients for the first q covariates.
//'@return a list of four variables
//'\describe{
//'  \item{x}{an n by p matrix, n observations of p-dimensional covariates}
//'  \item{y}{an n by 1 vector, n observations of response variable}
//'  \item{sigma2}{a positive number for noise variance}
//'  \item{beta_nonzero}{the true non-zero regression coefficients for the first q covariates}
//'}
//'@author Jian Kang \email{<jiankang@@umich.edu>}
//'@examples
//'dat <- fast_data_gen(1000,10000,0.9,c(1,2,3))
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
  double sigma2 = var_y*(1.0 - R2)/R2;
  y += arma::randn(n)*sqrt(sigma2);
  return Rcpp::List::create(Named("x") = x,
                            Named("y") = y,
                            Named("sigma2") = sigma2,
                            Named("beta_nonzero") = beta_nonzero);
}



//'Posterior mean variable screening based on local spatial smoothing
//'
//'@param x a design matrix of dimension n by p, where n is the sample size and p
//'is the number of spatially varying predictors
//'@param y a vector of length n for response variable
//'@param coords a matrix of dimension p by d (d>0) for the spatial locations,
//'where d is the dimension of space
//'@param neighbors a matrix of dimension p by q (q>0) for indices of neighbors
//'@param num_neighbors a vector of length p for the number of
//'neighbors for each predictor
//'@param rho: a vector of length p for spatial smoothness parameters
//'@param theta: positive number control the prior precision. Default value is 1e-3.
//'@return a list of two variables,
//'\describe{
//' \item{pms_select}{indices of the predictors from most important to least important}
//' \item{pms_stat}{a p by 1 vector of posterior mean statistics}
//'}
//'@author Jian Kang \email{<jiankang@@umich.edu>}
//'@examples
//' maskfile <- file.path(system.file("nifti", package="PMS"),"MNI-maxprob-thr0-2mm.nii.gz")
//' mask <- oro.nifti::readNIfTI(maskfile)
//' imgfiles <- file.path(system.file("nifti", package="PMS"),sprintf("VBM_example_0%d.nii.gz",1:5))
//' img_dat <- fast_read_imgs_mask(imgfiles,maskfile)
//' img1 <- oro.nifti::readNIfTI(imgfiles[1])
//' nb <- find_brain_image_neighbors(img1, mask, radius=1)
//' # simulate data
//' n = 500
//' p = ncol(img_dat)
//' x = matrix(rnorm(n*p),nrow=n,ncol=p)
//' beta_coef = rep(0,length=ncol(img_dat))
//' true_idx = nb$mask_img_nb[249,]
//' beta_coef[true_idx] = 1
//' y = x%*%beta_coef + rnorm(nrow(img_dat),sd=0.1)
//' nb_cor <- 0.99
//' rho <- rep(-log(nb_cor)/4,length=ncol(img_dat))
//' res <- fast_PMS_local_spatial(x=x, y = y, coords=nb$maskcoords,neighbors=nb$mask_img_nb,num_neighbors=nb$num_neighbors, rho = rho)
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
  arma::vec sqrt_rho = sqrt(rho);

  for(int j=0;j<x.n_cols;j++){
    for(int k=0;k<num_neighbors(j);k++){
        int l = neighbors(j,k) - 1L;
        arma::rowvec temp = (coords.row(l) - coords.row(j));
        double dist = arma::accu(temp%temp);
        lx.col(j) += exp(-sqrt_rho(j)*sqrt_rho(l)*dist)*x.col(l);
    }
  }

  arma::mat XXt = x*lx.t();
  XXt.diag() += theta;
  arma::vec screen_stat = arma::solve(XXt,y,arma::solve_opts::likely_sympd);
  screen_stat = lx.t()*screen_stat;
  arma::vec abs_screen_stat = arma::abs(screen_stat);
  arma::uvec od = arma::sort_index(abs_screen_stat,"descend");
  return Rcpp::List::create(Named("pms_select") = od+1L,
                            Named("pms_stat") = screen_stat);
}

//'Fast read multiple imaging data files on a brain mask
//'@param imgfiles a character vector of imaging file names
//'@param maskfile a character scalar for the  mask file
//'@param verbose a positive integer showing the progress on the number of images being loaded
//'@return a matrix of dimension n by p, where n is the number of images and p is the number of voxels on the mask
//'@author Jian Kang \email{<jiankang@@umich.edu>}
//'@examples
//' maskfile <- file.path(system.file("nifti", package="PMS"),"MNI-maxprob-thr0-2mm.nii.gz")
//' mask <- oro.nifti::readNIfTI(maskfile)
//' imgfiles <- file.path(system.file("nifti", package="PMS"),sprintf("VBM_example_0%d.nii.gz",1:5))
//' img_dat <- fast_read_imgs_mask(imgfiles,maskfile)
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

  return img_mat.t();
}




