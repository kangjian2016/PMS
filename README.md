# R package `PMS`


## Overview

PMS(Posterior Mean Screening) is a variable screening method constructed from the Bayesian framework. Different from other variable screening methods, the PMS method could incorporate both prior mean and prior covariance matrix information systematically. We construct the specific formulation of screening statistics under several different prior structures. 


## Install from Github

If the following R packages: `devtools`, `mcmc`, `truncnorm`, `glmnet` have not been installed, we should install them first:

```r
install.packages("devtools")
install.packages("mcmc")
install.packages("truncnorm")
install.packages("glmnet")
```
And then, install package `PMS` from Github with 

```r
devtools::install_github("kangjian2016/PMS")
library(PMS)
```
## Documentation

The code for documentation page:
```r
?pms_screening
?fast_PMS_cpp
?fast_PMS_local_spatial
?fast_read_imgs_mask
?find_brain_image_neighbors
```

## An example

```r
set.seed(1000)
rs = 9#signal to noise ratio
n = 200
p = 10000
theta = 1e-3
beta_nonzero = c(3,3,3,3,3,-7.5)
dat = data_gen(n,p,rs,beta_nonzero)
S_0 = c(2,6)
Lambda = Matrix::sparseMatrix(i=c(1:p),j=c(1:p),x = rep(1,p))
Lambda_s =  Matrix::sparseMatrix(i=c(1:2),j=c(1:2),x = rep(1,2))
PMS = pms_screening(x = dat$x,y = dat$y,family = "gaussian",
method = "selection",Lambda = Lambda,Lambda_s = Lambda_s,
theta = theta,selected_num = p,idx_s = S_0)$pms_select
```

## Another example for imaging data
```r
# read mask file
maskfile <- file.path(system.file("nifti", package="PMS"),"MNI-maxprob-thr0-2mm.nii.gz")
mask <- oro.nifti::readNIfTI(maskfile)
# read multiple image files on brain mask
imgfiles <- file.path(system.file("nifti", package="PMS"),sprintf("VBM_example_0%d.nii.gz",1:5))
img_dat <- fast_read_imgs_mask(imgfiles,maskfile)
#find neighboring voxels
img1 <- oro.nifti::readNIfTI(imgfiles[1])
nb <- find_brain_image_neighbors(img1, mask, radius=1)
# simulate data based on the image spatial structure
n = 500
p = ncol(img_dat)
x = matrix(rnorm(n*p),nrow=n,ncol=p)
beta_coef = rep(0,length=ncol(img_dat))
true_idx = nb$mask_img_nb[249,]
beta_coef[true_idx] = 1
y = x%*%beta_coef + rnorm(nrow(img_dat),sd=0.1)
nb_cor <- 0.99
rho <- rep(-log(nb_cor)/4,length=ncol(img_dat))
res <- fast_PMS_local_spatial(x=x, y = y, coords=nb$maskcoords,neighbors=nb$mask_img_nb,num_neighbors=nb$num_neighbors, rho = rho)

```
## Reference

Jie He and Jian Kang (2020) Prior knowledge guided ultra-high dimensional variable screening with application to neuroimaging data. Statistica Sinica, In Press. 

