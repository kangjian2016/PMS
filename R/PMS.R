#' Posterior Mean Variable Screening
#'
#' A variable screening method for linear regression model with different constructions of prior knowledge based on the framework of Bayesian modelling
#'
#'@export
#'@importFrom Matrix sparseMatrix
#'@importFrom Matrix bdiag
#'@importFrom truncnorm rtruncnorm
#'@importFrom mcmc metrop
#'@importFrom stats binomial
#'@importFrom stats coef
#'@importFrom stats gaussian
#'@importFrom stats poisson
#'@importFrom stats runif
#'@importFrom stats pnorm
#'@importFrom stats kmeans
#'@importFrom methods as
#'@importFrom glmnet glmnet
#'
#'@param x a numeric matrix represents the predictor variables with each row be a p-dimensional observation.
#'@param y a numeric vector of response with each component be an observation.
#'@param family the model type whose choices are the same as glmnet.
#'@param method the specific construction of prior knowlwdge in PMS method. "spatial_fixed" indicates Lambda is a known matrix;
#'              "selection", "rank" and "group" respectively stands for PMS method apply prior selected features set,
#'               ranking information as well as group structure as prior mean knowledge; and "spatial", "network" respectively means that prior
#'               covariance matrix Lambda has spatial and network structure with unknown parameters.
#'@param mu a numeric vector stands for known prior mean information, which can be simply set as "0" if no prior mean information is provided in methods "spatial" and "network".
#'@param Lambda a numeric matrix of known prior covariance knowledge. Similarly, we can choose "Lambda" as identity matrix if prior covariance knowledge is unavailable in
#'       "selection", "rank" or "group" methods.
#'@param Lambda_s a numeric matrix represents the correlation matrix of coefficients in prior selected features(different groups) when apply "selection" or "group" method.
#'@param theta a numeric scalar corresponds to the value of \eqn{\theta} in PMS statistics.
#'@param C_rho a numeric 2-dimensional vector with 2 components correspond to the endpoints of the intervel appears in the distribution of parameter rho in
#'             matrix Lambda when use "spatial" method.
#'@param L a numeric matrix stands for the normalized graph Laplacian matrix in Lambda under the "network" method.
#'@param B a numeric matrix stands for the group structure when apply "group" method.
#'@param idx_s a vecter corresponds to the index of prior selected features, which is only needed in "selection" and "rank" methods.
#'@param rank_s a vector corresponds to the ranking information of prior selected features in "rank" method.
#'@param xgrid a matrix represents the grids for covariates under image regression case.
#'@param selected_num an integer ranges from 1 to p, which represents the number of features selected.
#'@param num_1 an integer, describes the total times of sampling during parameter estimation in "rank", "spatial" and "network" methods.
#'@param cut_num_1 an integer, we throw away the first cut_num_1 samples when estimate parameter in "rank", "spatial" and "network" methods by sampling.
#'@param leng an integer stands for the sampling interval in the procedure of parameter estimation in "rank", "spatial" and "network" methods.
#'@param trho a numeric scalar controlling the degree of sparsity to matrix Lambda in method "spatial".
#'@return a list of three variables, where "pms_select" represents the index of selected features with number of selected_num,
#'        "pms_abs" is the absolute value of PMS statistics and "elapsed" be the computation time.
#'@references Prior knowledge guided ultra-high dimensional variable screening with application to neuroimaging data,working paper.
#'@author Jie He \email{<jiehe@@umich.edu>}, Jian Kang \email{<jiankang@@umich.edu>}
#'@examples
#'
#'set.seed(1000)
#'rs = 9
#'n = 200
#'p = 10000
#'theta = 1e-3
#'beta_nonzero = c(3,3,3,3,3,-7.5)
#'dat = data_gen(n,p,rs,beta_nonzero)
#'S_0 = c(2,6)
#'Lambda = Matrix::sparseMatrix(i=c(1:p),j=c(1:p),x=rep(1,p))
#'Lambda_s = Matrix::sparseMatrix(i=c(1:2),j=c(1:2),x=rep(1,2))
#'PMS_res = pms_screening(x = dat$x,y = dat$y,family = "gaussian",
#'method = "selection",Lambda = Lambda,Lambda_s = Lambda_s,
#'theta = theta,selected_num = p,idx_s = S_0)
#'
pms_screening = function (x, y, family = "gaussian", method = "spatial_fixed", mu = NULL ,Lambda = NULL,
                          Lambda_s = NULL, theta = NULL, C_rho = NULL, L = NULL, B = NULL, idx_s = NULL,
                          rank_s = NULL, xgrid = NULL, selected_num = NULL, num_1 = NULL, cut_num_1 = NULL,
                          leng = NULL, trho = NULL) {
  elapsed = proc.time()[3]
  X = as.matrix(x)
  Y = y
  if (is.null(dim(x))) {
    p = 1
    n = length(x)
  }else {
    p = dim(x)[2]
    n = dim(x)[1]
  }
  if (p == 1 || (p < n)) {
    selected_num = p
  }else {
    selected_num = selected_num
  }
  if(is.null(selected_num)){
    selected_num  = p
  }

  if(is.null(Lambda)){
    Lambda = Matrix::sparseMatrix(i=c(1:p),j=c(1:p),x=rep(1,p))
  }

  if(is.null(theta)){
    theta = 1e-3
  }

  if (family == "gaussian") {
     if (method == "spatial_fixed"){
       PMS_all = pms_spatial_fixed(x = X,y = Y,Lambda = Lambda, theta = theta)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else if (method == "selection"){
       PMS_all = pms_selection(x = X,y = Y,idx_s = idx_s, Lambda_s = Lambda_s, Lambda = Lambda, theta = theta)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else if (method == "rank"){
       PMS_all = pms_rank(x = X,y = Y,idx_s = idx_s, Lambda_s = Lambda_s, Lambda = Lambda, theta = theta,
                                     rank_s = rank_s, num_1 = num_1, cut_num_1 = cut_num_1, leng = leng)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else if (method == "group"){
       PMS_all = pms_group(x = X,y = Y,Lambda_s = Lambda_s, Lambda = Lambda, B = B,theta = theta)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else if (method == "spatial"){
       PMS_all = pms_spatial(x = X,y = Y, xgrid = xgrid, mu = mu, theta = theta, num_1 = num_1, cut_num_1 = cut_num_1, trho = trho, C_rho = C_rho)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else if (method == "network"){
       PMS_all = pms_network(x = X,y = Y,L = L, mu = mu, theta = theta, num_1 = num_1, cut_num_1 = cut_num_1)
       ranking = PMS_all$ranking
       pms_abs = PMS_all$PMS_0
     }else{
       stop("This method is unknown and not supported/")
     }
    }else {
      model = glmnet(x = X, y = Y, family = family, alpha = 0, lambda = 1, intercept = FALSE)
      coefs = coef(model)[-1]
      pms_abs = abs(coefs)
      ranking = order(pms_abs, decreasing = TRUE)
    }
    pms_select = ranking[1:selected_num]
  return(list(pms_select = pms_select, pms_abs = pms_abs, elapsed = proc.time()[3]-elapsed))#return feature index and time
}


inv_omega = function(x, theta = NULL, Lambda_s = NULL, B = NULL){
  n = dim(x)[1]
  p = dim(x)[2]
  x = x %*% Matrix::t(B)
  Lambda_x = x %*% Lambda_s
  inv_om = Lambda_x %*% Matrix::t(x)
  inv_om = as.matrix(inv_om) + theta * diag(n)
  return(inv_om)
}


mu_selection = function(x, y, idx_s = NULL, theta = NULL, Lambda_s = NULL){
  p = dim(x)[2]
  s = length(idx_s)
  B = matrix(0,s,p) #matrix to show which element to be selected
  for (i in 1:s){
    B[i,idx_s[i]] = 1
  }
  B = as(B, "sparseMatrix")
  inv_om = inv_omega(x = x, theta = theta, Lambda_s = Lambda_s, B = B)
  x_s = x[ ,idx_s]
  xy = cbind(x_s,y)
  om_xy = solve(inv_om, xy)
  om_x = om_xy[,1:s]
  om_y = om_xy[,(s+1)]
  K_1 = solve(t(x_s) %*% om_x + 1e-3 * diag(s))
  mu_s = K_1 %*% t(x_s) %*% om_y
  mu = rep(0, p)
  mu[idx_s] = mu_s
  return(list(mu = mu, K = K_1))
}


mean_sd_prob = function(index, mu_s, nu_s, K_s, order_s){
  q = length(nu_s)
  K_0 = K_s[-order_s[index],-order_s[index]]
  dim_0 = ncol(K_0)
  if (is.null(dim_0)==1){
    if (K_0 < 1e-03){
      K_0 = K_0 + 1e-03
    }
    M_0 = 1/K_0
  }else{
    M_0 = solve(K_0 + 1e-3 * diag(q-1))
  }
  M_1 = K_s[order_s[index], -order_s[index]]
  K_1 = M_1 %*% M_0 %*% (mu_s[-order_s[index]] - nu_s[-order_s[index]])
  K_2 = K_s[order_s[index], order_s[index]] - t(M_1) %*% M_0 %*% M_1
  mean_1 = nu_s[order_s[index]] - K_1
  sd_1 = sqrt(abs(K_2 + 1e-06))
  if (index > 1 & index < q){
    log_w1 = log(pnorm(abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) - pnorm(abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) + 1e-6)
    log_w2 = log(pnorm(abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) - pnorm(abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) + pnorm(-abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) - pnorm(-abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) + 1e-6)
  } else if (index == 1) {
    log_w1 = log(1 - pnorm(abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) + 1e-6)
    log_w2 = log(1 - pnorm(abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) + pnorm(-abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1) + 1e-6)
  } else {
    log_w1 = log(pnorm(abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) - pnorm(0, mean = mean_1, sd = sd_1) + 1e-6)
    log_w2 = log(pnorm(abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) - pnorm(-abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1) + 1e-6)
  }
  return(list(mean_1 = mean_1, sd_1 = sd_1, prob = exp(log_w1 - log_w2)))
}


mu_update = function(index, mu_s, nu_s, K_s, order_s){
  q = length(nu_s)
  m_s_p = mean_sd_prob(index = index, mu_s = mu_s, nu_s = nu_s, K_s = K_s, order_s = order_s)
  mean_1 = m_s_p$mean_1
  sd_1 = m_s_p$sd_1
  prob_1 = m_s_p$prob
  prob_1 = round(prob_1,2)
  indicator = sample(1:2, 1, prob = c(prob_1, 1 - prob_1))
  if (indicator == 1){
    if (index > 1 & index < q){
      mu_new = rtruncnorm(1, a = abs(mu_s[order_s[index+1]]), b = abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1)
    } else if(index == 1){
      mu_new = rtruncnorm(1, a = abs(mu_s[order_s[index+1]]), b = 1e03, mean = mean_1, sd = sd_1)#need package "truncnorm"
    } else {
      mu_new = rtruncnorm(1, a = 0, b = abs(mu_s[order_s[index-1]]), mean = mean_1, sd = sd_1)
    }
  } else{
    if (index > 1 & index < q){
      mu_new = rtruncnorm(1, a = - abs(mu_s[order_s[index-1]]),  b = - abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1)
    } else if(index == 1){
      mu_new = rtruncnorm(1, a = -1e03, b = - abs(mu_s[order_s[index+1]]), mean = mean_1, sd = sd_1)
    } else {
      mu_new = rtruncnorm(1, a = - abs(mu_s[order_s[index-1]]), b = 0, mean = mean_1, sd = sd_1)
    }
  }
  return(mu_new = mu_new)
}


mu_gibbs = function(nu_s, K_s, num, cut_num, leng, order_s){
  q = length(nu_s)
  mu_sample = NULL
  mu_0 = runif(q,-10,10)
  order_0 = order(abs(mu_0), decreasing = TRUE)
  mu_0[order_s] = mu_0[order_0]
  if (q >= 3){
    for(i in 1:num){
      mu_0[order_s[1]] = mu_update(index = 1, mu_s = mu_0, nu_s = nu_s, K_s = K_s, order_s = order_s)
      for (j in 2:(q-1)){
        mu_1 = mu_0[order_s[j-1]]
        mu_2 = mu_0[order_s[j+1]]
        if (abs(mu_1 - mu_2) < 1e-6){
          mu_n = 1/2* (mu_1 + mu_2)
        } else{
          mu_n = mu_update(index = j, mu_s = mu_0, nu_s = nu_s, K_s = K_s, order_s = order_s)
        }
        mu_0[order_s[j]] = mu_n
      }
      mu_0[order_s[q]] = mu_update(index = q, mu_s = mu_0, nu_s = nu_s, K_s = K_s, order_s = order_s)
      if ( (i > cut_num)  &  ((i- cut_num) %% leng == 0)){
        mu_sample = cbind(mu_sample, mu_0)
      }
    }
  } else if (q<=2){
    for(i in 1:num){
      mu_1 = mu_update(index = 1, mu_s = mu_0, nu_s = nu_s, K_s = K_s, order_s = order_s)
      mu_0[order_s[1]] = mu_1
      mu_2 = mu_0[order_s[2]]
      if (abs(mu_1-mu_2) < 1e-6){
        mu_n = 1/2 *(mu_1 + mu_2)
      } else{
        mu_n = mu_update(index = 2, mu_s = mu_0, nu_s = nu_s, K_s = K_s, order_s = order_s)
      }
      mu_0[order_s[2]] = mu_n
      if ( (i > cut_num)  &  ((i- cut_num) %% leng == 0)){
        mu_sample = cbind(mu_sample, mu_0)
      }
    }
  }
  return(mu_sample)#a q by num matrix
}


mu_ranking = function(x, y, idx_s = NULL, theta = NULL, Lambda_s = NULL, rank_s = NULL, num = NULL, cut_num = NULL, leng = NULL){
  p = dim(x)[2]
  order_s = order(rank_s, decreasing = TRUE)#the index from important to unimportant
  nu_K = mu_selection(x = x, y = y, idx_s = idx_s, theta = theta, Lambda_s = Lambda_s)
  nu_s = nu_K$mu[idx_s]
  K_s = nu_K$K
  mu_hat_mat = mu_gibbs(nu_s = nu_s, K_s = K_s, num = num, cut_num = cut_num,
                        leng = leng, order_s = order_s)
  mu_hat_0 = rowMeans(mu_hat_mat)
  mu_hat = rep(0, p)
  mu_hat[idx_s] = mu_hat_0
  return(mu_hat)
}


mu_group = function(x, y, B = NULL, Lambda_s = NULL, theta = NULL){
  B = as(B, "sparseMatrix")
  inv_om = inv_omega(x = x, Lambda_s = Lambda_s, theta = theta, B = B)
  x_b = x %*% Matrix::t(B)
  xy = cbind(x_b, y)
  om_xy = solve(inv_om, xy)
  om_x = om_xy[,1:ncol(x_b)]
  om_y = om_xy[,(ncol(x_b)+1)]
  M_1 = t(x_b) %*% om_x + (1e-3) *diag(nrow(B))
  mu_b = solve(M_1, t(x_b) %*% om_y)
  mu_group = Matrix::t(B) %*% mu_b
  mu_group = as.vector(mu_group)
  return(mu_group)
}



rho_posterior = function(rho, xgrid = NULL, x = NULL , y = NULL, trho = NULL, theta = NULL, mu = NULL, C_rho = NULL){
  p = ncol(x)
  Lambda = Lambda_sparse(xgrid = xgrid, rho = rho, trho = trho, num_group = 10)
  B_1 = sparseMatrix(i= c(1:p), j = c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x, theta = theta, Lambda_s = Lambda, B = B_1)
  inv_om = 1/2*(inv_om + Matrix::t(inv_om))
  Omega = solve(inv_om)
  om_val_0 = eigen(inv_om)$values
  if (min(om_val_0)< 0){
    om_val_0 = om_val_0 + abs(min(om_val_0)) + 1e-3
  }
  om_val = om_val_0^(-1)
  log_om_det = (1/2) * sum(log(om_val))
  res = y - x %*% mu
  ind_rho = (rho <= C_rho[2] & rho > C_rho[1])
  if (ind_rho == 1){
    posterior_0 = (1/2)*t(res) %*% Omega %*% res
    if (abs(posterior_0) <= 1e04){
      posterior_rho = log_om_det - (1/2)*posterior_0
    }else{
      posterior_rho = -1e04
    }
  }else{
    posterior_rho = -1e04
  }
  return(posterior_rho)
}



Lambda_spatial = function(x, y, xgrid, mu = NULL, C_rho = NULL, theta = NULL, num = NULL, cut_num = NULL, trho = NULL){
  rho_init = runif(1,C_rho[1],C_rho[2])
  rho_sample = metrop(obj = rho_posterior, initial = rho_init, nbatch = num, xgrid = xgrid, trho = trho, theta = theta, mu = mu, x = x, y = y, C_rho = C_rho)$batch
  if (cut_num == 0){
    rho_sample = rho_sample
  }else{
    idx_0 = c(1:cut_num)#the first cut_num samples to be thorw away
    rho_sample = rho_sample[-idx_0]
  }
  k = length(rho_sample)
  Lambda_list = lapply(1:k, function(i) Lambda_sparse(xgrid = xgrid, rho = rho_sample[i], trho = trho, num_group = 10))
  Lambda_spatial = (1/k) * Reduce('+', Lambda_list)
  #Lambda_spatial = as.matrix(Lambda_spatial)
  rm(Lambda_list)
  return(Lambda_spatial)
}



Lambda_network_0 = function(L, epsi){
  p = ncol(L)
  I_0 = sparseMatrix(i = c(1:p), j = c(1:p), x = rep(1,p))
  L_1 = L + epsi * I_0
  Lambda_network_0 = solve(L_1)
  return(Lambda_network_0)
}


epsi_posterior = function(epsi, x = NULL, y = NULL, theta = NULL, L = NULL, mu = NULL){
  p = dim(x)[2]
  Lambda = Lambda_network_0(L = L, epsi = epsi)
  B_1 = sparseMatrix(i = c(1:p), j = c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x,theta = theta, Lambda_s = Lambda, B = B_1)
  inv_om = 1/2 * (inv_om + t(inv_om))
  Omega = solve(inv_om)
  om_val_0 = eigen(inv_om)$values
  if (min(om_val_0) < 0){
    om_val_0 = abs(min(om_val_0)) + 1e-3
  }
  om_val = om_val_0^(-1)
  log_om_val = sum(log(om_val))
  res = y - x %*% mu
  if(epsi > 0){
    posterior_0 = t(res) %*% Omega %*% res
    #posterior_0 = as.numeric(posterior_0)
    if (abs(posterior_0) <= 1e04){
      posterior_epsi = log_om_val - (1/2) * posterior_0
    }else{
      posterior_epsi = -1e04
    }
  }else{
    posterior_epsi = -1e04
  }
  return(posterior_epsi)
}


Lambda_network = function(x, y, mu = NULL, L = NULL, theta = NULL, num = NULL,
                          cut_num = NULL){
  epsi_0 = runif(1, -3, -2)
  epsi_init = 10^epsi_0
  epsi_sample = metrop(obj = epsi_posterior, initial = epsi_init, nbatch = num, x = x, y = y, theta = theta, L = L, mu = mu)$batch
  if(cut_num == 0){
    epsi_sample = epsi_sample
  }else{
    idx_0 = c(1:cut_num)
    epsi_sample = epsi_sample[-idx_0]
  }
  k = length(epsi_sample)
  Lambda_list = lapply(1:k, function(i) Lambda_network_0 (L = L, epsi = epsi_sample[i]))
  Lambda_network = (1/k) * Reduce("+", Lambda_list)
  return(Lambda_network)
}


Lambda_sparse = function(xgrid, rho = NULL, trho = NULL, num_group = NULL){
  d = ncol(xgrid)
  group_a = kmeans(x = xgrid, centers = num_group)$cluster
  n_group = length(unique(group_a))
  group_all = NULL
  for (i in 1:n_group){
    group_list = which(group_a == i)
    p_n = length(group_list)
    group_all = c(group_all, group_list)
    index_c = NULL
    index_r = NULL
    for (j in 1:p_n){
      index_a = c(1:j)
      index_b = rep(j,j)
      index_c = c(index_c,index_a)
      index_r = c(index_r,index_b)
    }
    distant_0 = lapply(1:d, function(i) (xgrid[group_list[index_c],i]-xgrid[group_list[index_r],i])^2)
    distant_1 = Reduce("+", distant_0)
    index_final = which(distant_1 <= -(log(trho)/rho))
    index_r_final = index_r[index_final]
    index_c_final = index_c[index_final]
    Elements = exp(-rho*distant_1[index_final])
    Lambda_0 = sparseMatrix(i=index_r_final, j = index_c_final, x = Elements)
    Lambda_1 = t(Lambda_0)
    diag(Lambda_1) = 0
    Lambda_0 = Lambda_0 + Lambda_1
    if (i==1){
      Lambda = Lambda_0
    }else{
      Lambda = bdiag(Lambda,Lambda_0)
    }
  }
  group_order = order(group_all)
  Lambda = Lambda[group_order, group_order]
  return(Lambda)
}


pms_spatial_fixed = function(x,y,Lambda = NULL, theta = NULL){
  p = ncol(x)
  mu = rep(0,p)
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  Mat = solve(inv_om, y)
  PMS = Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0,decreasing = TRUE)
  return(list(ranking = ranking,PMS_0 = PMS_0))
}

pms_selection = function(x,y, idx_s = NULL, Lambda_s = NULL, Lambda = NULL, theta = NULL){
  p = ncol(x)
  mu = mu_selection (x, y, idx_s = idx_s, theta = theta, Lambda_s = Lambda_s)$mu
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  res = y - x %*% mu
  Mat = solve(inv_om, res)
  PMS = mu + Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0, decreasing = TRUE)
  return(list(ranking = ranking, PMS_0 = PMS_0))
}

pms_rank = function(x,y, idx_s = NULL, Lambda_s = NULL, Lambda = NULL, theta = NULL,
                    rank_s = NULL, num_1 = NULL, cut_num_1 = NULL, leng = NULL){
  p = ncol(x)
  mu = mu_ranking(x, y, idx_s = idx_s, theta = theta, Lambda_s = Lambda_s, rank_s = rank_s, num = num_1, cut_num = cut_num_1, leng = leng)
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  res = y - x %*% mu
  Mat = solve(inv_om, res)
  PMS = mu + Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0, decreasing = TRUE)
  return(list(ranking = ranking, PMS_0 = PMS_0))
}

pms_group = function(x,y, Lambda_s = NULL, Lambda = NULL, B = NULL,theta = NULL){
  p = ncol(x)
  mu = mu_group(x, y,B = B, Lambda_s = Lambda_s, theta = theta)
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  res = y - x %*% mu
  Mat = solve(inv_om, res)
  PMS = mu + Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0, decreasing = TRUE)
  return(list(ranking = ranking, PMS_0 = PMS_0))
}

pms_spatial = function(x,y, xgrid = NULL, mu = NULL, theta = NULL, num_1 = NULL, cut_num_1 = NULL, trho = NULL, C_rho = NULL){
  p = ncol(x)
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  Lambda = Lambda_spatial(x,y,xgrid = xgrid, mu = mu, C_rho = C_rho, theta = theta, num = num_1, cut_num = cut_num_1, trho = trho)
  #Lambda = as(Lambda,"sparseMatrix")
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  res = y - x %*% mu
  Mat = solve(inv_om, res)
  PMS = mu + Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0, decreasing = TRUE)
  return(list(ranking = ranking, PMS_0 = PMS_0))
}

pms_network = function(x,y,L = NULL, mu = NULL, theta = NULL, num_1 = NULL, cut_num_1 = NULL){
  p = ncol(x)
  B_1 = sparseMatrix(i = c(1:p), j =c(1:p), x = rep(1,p))
  Lambda = Lambda_network(x,y,mu = mu,L = L, theta = theta, num = num_1, cut_num = cut_num_1)
  inv_om = inv_omega(x = x, Lambda_s = Lambda, theta = theta, B = B_1)
  x_n = x %*% Lambda
  res = y - x %*% mu
  Mat = solve(inv_om, res)
  PMS = mu + Matrix::t(x_n) %*% Mat
  PMS_0 = abs(PMS)
  ranking = order(PMS_0, decreasing = TRUE)
  return(list(ranking = ranking, PMS_0 = PMS_0))
}





#' Find the neighbors for each voxel in images
#' @param img an nifti object.
#' @param mask an nifti object.
#' \code{mask > 0} specifies which voxels are on the mask.
#' Default value is NULL indicating all voxels are considered.
#' @param radius the size (voxel) of neighborhood.
#' @return the neighbor indices of each voxel in each row.
#' @author Jian Kang <jiankang@umich.edu>
#' @examples
#' maskfile <- file.path(system.file("nifti", package="PMS"),"AAL_MNI_2mm.nii")
#' mask <- oro.nifti::readNIfTI(maskfile)
#' imgfile <- file.path(system.file("nifti", package="PMS"),"VBM_example.nii.gz")
#' img <- oro.nifti::readNIfTI(imgfile)
#' nb <- find_brain_image_neighbors(img, mask,radius=1)
#' @export
#'
find_brain_image_neighbors <- function(img, mask=NULL, radius = 1){


  grids <- list(X = 1:img@dim_[2],
                       Y = 1:img@dim_[3],
                       Z = 1:img@dim_[4])

  d = length(grids)
  dim = sapply(1:d,function(i) length(grids[[i]]))
  coords = expand.grid(grids)
  if(!is.null(mask)){
    maskidx <- which(mask>0)
    num_voxels = length(maskidx)
  } else{
    num_voxels = prod(dim)
  }


  nb_idx = seq(-radius,radius,by=1)
  nb_idx_list = list()
  for(i in 1:d){
    nb_idx_list[[i]] = nb_idx
  }
  idx_patterns = expand.grid(nb_idx_list)
  zero_idx = which(apply(idx_patterns==0,1,all))
  idx_patterns = idx_patterns[-zero_idx,]
  nb = array(NA,dim=c(num_voxels,(2*radius+1)^d))
  #nb_dist = array(0, dim=c(num_voxels,(2*radius+1)^d))
  if(!is.null(mask)){
    nb[,1] = maskidx
  } else{
    nb[,1] = 1:num_voxels
  }
  pb = txtProgressBar(style=3)
  for(l in 1:nrow(idx_patterns)){
    img_arr_idx = list()
    nb_arr_idx = list()
    for(i in 1:d){
      img_arr_idx[[i]] = max(1-idx_patterns[l,i],1):min(dim[i] - idx_patterns[l,i],dim[i])
      nb_arr_idx[[i]] = max(1+idx_patterns[l,i],1):min(dim[i] + idx_patterns[l,i],dim[i])
    }
    img_arr = expand.grid(img_arr_idx)
    img_vec_idx = 0
    for(i in d:1){
      img_vec_idx = dim[i]*img_vec_idx + (img_arr[,i]-1)
    }
    img_vec_idx = img_vec_idx + 1

    nb_arr = expand.grid(nb_arr_idx)
    nb_vec_idx = 0
    for(i in d:1){
      nb_vec_idx = dim[i]*nb_vec_idx + (nb_arr[,i]-1)
    }
    nb_vec_idx = nb_vec_idx + 1
    if(!is.null(mask)){
      nb_vec_idx_0 = intersect(nb_vec_idx,maskidx)
      subidx = match(nb_vec_idx_0,nb_vec_idx)
      nb_vec_idx = match(nb_vec_idx_0,maskidx)
      img_vec_idx = img_vec_idx[subidx]
    }
    nb[nb_vec_idx,l+1] = img_vec_idx
    #for(s in 1:d){
    #  nb_dist[,l+1] = nb_dist[,l+1] + (coords[nb[,l+1],s] - coords[,s])^2
    #}
    setTxtProgressBar(pb,l/nrow(idx_patterns))
  }

  if(!is.null(mask)){
    for(l in 2:ncol(nb)){
      nb[,l] = ifelse(is.element(nb[,l],maskidx),nb[,l],NA)
    }
  }

  close(pb)
  return(nb)
}

