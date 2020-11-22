# R package `PMS`


## Overview

PMS(Posterior Mean Screening) is a variable screening method constructed from the Bayesian framework. Different from other variable screening methods, the PMS method could incorporate both prior mean and prior covariance matrix information systematically. We construct the specific formulation of screening statistics under several different prior structures. 


## Install from Github

If R package `devtools` has not been installed, we should install it firstly:

```r
install.packages("devtools")
```
And then, install package `PMS` from Github with 

```r
devtools::install_github("https://github.com/kangjian2016/PMS.git")
libaray(PMS)
```
## Documentation

The code for documentation page:
```r
?pms_screening
```

## An example

```r
set.seed(1000)
rs = 9#signal to noise ratio
n=100
p=10000
theta = 1e-3
beta_nonzero = c(3,3,3,3,3,-7.5)
dat = data_gen(n,p,rs,beta_nonzero)
S_0 = c(2,6)
Lambda = Matrix::sparseMatrix(i=c(1:p),j=c(1:p),x = rep(1,p))
Lambda_s =  Matrix::sparseMatrix(i=c(1:2),j=c(1:2),x = rep(1,2))
PMS=pms_screening(x = dat$x,y = dat$y,family = "gaussian",method = "selection",Lambda = Lambda,Lambda_s = Lambda_s,
          theta = theta,selected_num = p,idx_s = S_0)$pms_select
```
## Reference

Jie He and Jian Kang (2020+), Prior Knowledge Guided Ultra-high Dimensional Variable Screening with Application to Neuroimaging Data. 

