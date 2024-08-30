
#### Libraries ####
library(splines)
library(MASS)
library(fda)
library(ddalpha)
library(emdbook)
library(RcppArmadillo)
library(splines2)
library(mvnfast)
library(coda)
library(LearnBayes)
library(lattice)

#### Rcpp functions ####
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends="RcppArmadillo")
Rcpp::cppFunction("arma::vec Arma_mvrnorm(const arma::vec& mean, const arma::mat& cov) {
  return arma::mvnrnd(mean, cov);}", depends="RcppArmadillo")
Rcpp::cppFunction("double  Arma_rnorm(const double &mu, const double &sigma){
      return arma::randn(arma::distr_param(mu,sigma));} ", depends="RcppArmadillo") #variance = sigma^2
Rcpp::cppFunction("double  Arma_runif(const double &a, const double &b){
      return arma::randu(arma::distr_param(a,b));} ", depends="RcppArmadillo")


#####

b_splines = function(n_knots, time_vec, coeff_vec=NULL,cond_order=F){
  if(cond_order){order = 4}else{order = 3}
  n_knots = n_knots-order
  knots = seq(0, 1, length.out = n_knots+2)[-c(1,n_knots+2)]
  B = bs(time_vec, knots = knots, intercept=cond_order)
  if(!is.null(coeff_vec)){
    return(B%*%coeff_vec)
  }
  else{
    return(B)
  }
}

omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  return (K)
}

stack_Matrix<-function(x_list,a=NULL){
  n=length(x_list)
  X=NULL
  if(!is.null(a)){
    for (i in 1:n){
      X= rbind(X,x_list[[i]]*a[i])
    }
  }
  else{
    for (i in 1:n){
      X= rbind(X,x_list[[i]])
    }
  }
  return(X)
}

pi_phi_Tel = function (phi, beta, gamma_i,sigma_eps, P, sigma_phi, Upsilon, Bh, y_i,cond_order_beta=F,cond_order_gamma=F){
  warped_t = Bh%*%phi
  prod1 = b_splines(length(beta), warped_t, beta,cond_order = cond_order_beta)
  prod2 = b_splines(length(gamma_i), warped_t, gamma_i,cond_order = cond_order_gamma)
  add1 = y_i - ( prod1 + prod2 )
  term1 = crossprod(add1)/sigma_eps
  term2 = t(phi-Upsilon)%*%P%*%(phi-Upsilon)/sigma_phi
  return(-0.5*(term1+term2))
}

pi_phi_GN = function (phi, beta, gamma_i,sigma_eps,csi,r,f, Bh, y_i,cond_order_beta = F,cond_order_gamma=F){
  warped_t = Bh%*%phi
  prod1 = b_splines(length(beta), warped_t, beta, cond_order_beta)
  prod2 = b_splines(length(gamma_i), warped_t, gamma_i, cond_order_gamma)
  add1 = y_i - ( prod1 + prod2 )
  term1 = as.vector(crossprod(add1)/sigma_eps)
  term2 = as.vector(dgamma(csi,shape = r,rate = f))
  return(-0.5*term1+log(term2))
}


interp_spline <- function(x, y, nout = length(y)) {
  ind_out <- seq(min(x), max(x), len = nout)
  spfit <- splinefun(x, y)
  return(spfit(ind_out))
}

JARA_smooting_groups = function(y,n_groups = 3, n_per_group = c(32,32,36), nburn=100,niter=100,n_knots_m =20,n_knots_gamma=10,a_eps = 3,
                                b_eps = .1,a_lambda = 3, b_lambda = .1,
                                a_gamma=3,b_gamma=0.1, cond_order_beta = F, cond_order_gamma=T, pre_smoot = F){
  
  n_obs_max  = dim(y)[1] 
  n_patients = dim(y)[2] # y has patients in the columns
  
  y_list = list()     # list of lists of groups of lists of individulas for group
  y_val  = list()     # list of lists of numerical values matrices for group 
  n_obs  = list()
  
  n0 = 1 
  nF = n_per_group[1]
  
  for(g in 1:n_groups){
    
    y_temp = y[,n0:nF]
    y_list_temp = lapply(seq_len(ncol(y_temp)), function(i) na.omit(y_temp[, i]))
    n_obs_temp=sapply(y_list_temp,length)
    y_temp = unlist(y_list_temp)
    y_temp = as.numeric(y_temp)
    
    y_list[[g]] = y_list_temp 
    y_val[[g]]  = y_temp 
    n_obs[[g]]  = n_obs_temp
    
    n0 = nF + 1 
    nF = nF + n_per_group[g+1]
  }
  
  rm(n0)
  rm(nF)
  
  if(cond_order_beta){order_beta = 4}else{order_beta = 3}
  if(cond_order_gamma){order_gamma = 4}else{order_gamma = 3}
  
  # equal to 4 for cubic splines
  p = n_knots_m + order_beta
  k = n_knots_gamma + order_gamma
  Omega = omega_P(p)
  
  time_group = list()
  
  for(g in 1:n_groups){
    time_group[[g]] = lapply(1:n_per_group[g], function(i)  seq(0,1, length.out = n_obs[[g]][i]))
  }
  
  # knots and Bspline matrix for common shape function m
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  knots_gamma = seq(0, 1, length.out = n_knots_gamma+2)[-c(1,n_knots_gamma+2)]
  
  Bm_beta_group  = list()
  Bm_gamma_group = list()
  
  for(g in 1:n_groups){
    Bm_beta_group[[g]] <- lapply(1:n_per_group[g], function(i) bs(time_group[[g]][[i]], knots = knots_m, intercept = cond_order_beta))
    Bm_gamma_group[[g]] <- lapply(1:n_per_group[g], function(i) bs(time_group[[g]][[i]], knots = knots_gamma, intercept = cond_order_gamma))
  }
  
  beta_0 <- rep(0,p)
  gamma_0<- rep(0,k)
  
  
  #----- Save Structures -----#
  
  nrun <- nburn + niter
  
  #-- Sigma epsilon --#
  sigma_eps_save <- numeric(nrun)
  sigma_eps_save[1] <- 1/rgamma(1,a_eps,rate= b_eps) 
  
  #-- Lambda --#
  lambda_save <- numeric(nrun)
  lambda_save[1] <-1/rgamma(1,a_lambda,rate= b_lambda) 
  
  #-- Beta --#
  beta_save_group <- list()
  for(g in 1:n_groups){
    beta_save_group[[g]] = matrix(NA, nrow = nrun, ncol = p)
    beta_save_group[[g]][1,] <- mvrnorm(n = 1, mu = beta_0, Sigma = solve(Omega/lambda_save[1]) )
  }
  
  X <- lapply(1:n_groups, function(g) stack_Matrix(Bm_beta_group[[g]]))
  
  #-- Sigma gamma --#
  sigma_gamma_save <- numeric(nrun)
  sigma_gamma_save[1] <- 1/rgamma(1,a_gamma,rate= b_gamma) 
  
  #-- Gamma --#
  gamma_save_group=list() 
  
  sigma_gamma_matrix=diag(k)*sigma_gamma_save[1]
  inv_sigma_gamma_matrix = solve(sigma_gamma_matrix)
  
  for(g in 1:n_groups){
    gamma_save_gruppo_G = list()
    for( i in 1:n_per_group[g]){
      gamma_save_ind = matrix(NA,nrun,k)
      gamma_save_ind[1,] = mvrnorm(n = 1, mu = gamma_0, Sigma = sigma_gamma_matrix )
      gamma_save_gruppo_G[[i]] = gamma_save_ind
    }
    gamma_save_group[[g]] = gamma_save_gruppo_G
  }
  
  rm(gamma_save_gruppo_G)
  rm(gamma_save_ind)
  
  #----- MCMC Loop -----#
  
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(iter in 2:nrun){
    
    setTxtProgressBar(bar, iter)
    
    #-- Update Beta --#
    
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = lapply(1:n_groups, function(g) inv_sigma_beta + crossprod(X[[g]])/(sigma_eps_save[iter-1]) )    
    
    V_beta = lapply(1:n_groups, function(g)   armaInv(inv_V_beta[[g]])  )
    
    vec_list= lapply(1:n_groups, function(g)  lapply(1:n_per_group[g], function(i) gamma_save_group[[g]][[i]][iter-1,]) )    
    G = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { matrice %*% vettore}, Bm_gamma_group[[g]], vec_list[[g]]))   
    C_gamma<- lapply(1:n_groups, function(g) rep(do.call(rbind, G[[g]]),1)) 
    
    m_beta = lapply(1:n_groups, function(g)   V_beta[[g]]%*%(t(X[[g]]) %*% ((y_val[[g]]-C_gamma[[g]])/sigma_eps_save[iter-1]) + inv_sigma_beta %*% beta_0 )  ) 
    Beta_udpate =  lapply(1:n_groups, function(g) Arma_mvrnorm(m_beta[[g]], V_beta[[g]]) ) 
    
    aggiorna_matrice <- function(matrice, vettore) {
      matrice[iter, ] <- vettore
      return(matrice)
    }
    
    beta_save_group <- mapply(aggiorna_matrice,  beta_save_group,  Beta_udpate , SIMPLIFY = FALSE)
    
    #-- Update lambda --#
    
    a_star = 0.5*p*n_groups + a_lambda
    beta_somma = lapply(1:n_groups, function(g) t(beta_save_group[[g]][iter,])%*% Omega %*%beta_save_group[[g]][iter,] )
    sum_term = as.numeric( sum(unlist(beta_somma)) )
    b_star = 0.5 * sum_term +  b_lambda
    lambda_save[iter]<-1/rgamma(n=1, shape=a_star, rate=b_star)
    
    #-- Update Gamma --#
    
    prod_X_tilde =  lapply(1:n_groups, function(g) lapply(Bm_gamma_group[[g]],function(matrice){crossprod(matrice)/sigma_eps_save[iter-1]})   )   
    inv_V_gamma  =  lapply(1:n_groups, function(g) lapply(prod_X_tilde[[g]],function(matrice){matrice + inv_sigma_gamma_matrix})  ) 
    V_gamma      =  lapply(1:n_groups, function(g) lapply(inv_V_gamma[[g]], function(matrice){armaInv(matrice)}) ) 
    
    c_tilde   = lapply(1:n_groups, function(g) lapply(Bm_beta_group[[g]],function(matrice){matrice%*%beta_save_group[[g]][iter,]})  ) 
    diff_term = lapply(1:n_groups, function(g) Map('-',y_list[[g]],c_tilde[[g]]) )
    prod_term = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { t(matrice) %*% vettore}, Bm_gamma_group[[g]], diff_term[[g]]) )
    m_gamma   = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { (1/sigma_eps_save[iter-1])*matrice %*% vettore}, V_gamma[[g]], prod_term[[g]]) ) 
    
    G = lapply(1:n_groups, function(g) Map(function(matrice, vettore) {mvrnorm(n = 1, mu = vettore, Sigma = matrice )},
                                           V_gamma[[g]], m_gamma[[g]]) ) 
    
    gamma_save_group <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g],
                                                              function(i) aggiorna_matrice(gamma_save_group[[g]][[i]],  G[[g]][[i]] ) ) ) 
    
    #-- Update Sigma gamma --#
    
    a_star=0.5*n_patients*k+a_gamma
    
    vettore_somma = lapply(1:n_groups, function(g)  lapply(1:n_per_group[g],function(i){crossprod(gamma_save_group[[g]][[i]][iter,])}) )
    b_star = b_gamma + 0.5 * sum(unlist(vettore_somma))
    
    sigma_gamma_save[iter]=1/rgamma(1,a_star,rate=b_star)
    
    inv_sigma_gamma_matrix=diag(k)/sigma_gamma_save[iter]
    
    
    #-- Update sigma eps --#
    
    a_star = a_eps + 0.5*sum(unlist(n_obs))
    
    term_beta = lapply(1:n_groups, function(g) lapply(Bm_beta_group[[g]],function(bm){bm%*%beta_save_group[[g]][iter,]}) )
    
    #vec_list è gia calcolato da prima 
    #vec_list= lapply(1:n_groups, function(g)  lapply(1:n_per_group[g], function(i) gamma_save_group[[g]][[i]][iter-1,]) ) 
    term_gamma = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { matrice %*% vettore}, Bm_gamma_group[[g]], vec_list[[g]]) )
    
    m_i = lapply(1:n_groups, function(g) Map("+", term_beta[[g]], term_gamma[[g]]) )
    diff_term = lapply(1:n_groups, function(g) Map("-", y_list[[g]], m_i[[g]]) ) 
    prod_term = lapply(1:n_groups, function(g) lapply(diff_term[[g]],function(v){crossprod(v)}) )
    somma=0.5*sum(unlist(prod_term))
    b_star=somma+b_eps
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
  }
  close(bar)
  
  beta_post_group  = lapply( 1: n_groups, function(g) colMeans(beta_save_group[[g]][-c(1:nburn),]) )
  
  gamma_post_group = lapply(1: n_groups, function(g) lapply(1:n_per_group[g], 
                                                            function(i) colMeans(gamma_save_group[[g]][[i]][-c(1:nburn),]) ))
  
  lambda_post <- mean(lambda_save[-c(1:nburn)])
  sigma_eps_post<-mean(sigma_eps_save[-c(1:nburn)])
  sigma_gamma_post<-mean(sigma_gamma_save[-c(1:nburn)])
  
  Bm_post_beta_group  = Bm_beta_group 
  Bm_post_gamma_group = Bm_gamma_group
  
  y_p <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) 
    as.numeric(Bm_post_beta_group[[g]][[i]] %*% beta_post_group[[g]] + Bm_post_gamma_group[[g]][[i]] %*% gamma_post_group[[g]][[i]])) ) 
  
  if(pre_smoot){
    print('PRE SMOOTING DONE')
    return(y_p)
  } else{
    print('DONE') 
    return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                            sigma_eps = sigma_eps_post, sigma_gamma=sigma_gamma_post,
                            gamma_group=gamma_post_group, y=y_p), 
                  full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                            gamma_group_s=gamma_save_group, sigma_gamma_s=sigma_gamma_save),
                  n_obs = n_obs,
                  Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group)))
  }
  
}




JARA_warping_GROUP = function(y,n_groups = 3, n_per_group = c(32,32,36), nburn=100,niter=100,
                              n_knots_m =20,n_knots_gamma=7,n_knots_h = 10,
                              a_eps = 3, b_eps = .1,a_lambda = 3, b_lambda = .1,a_lambda_gamma = 10, b_lambda_gamma = 2,
                              a_gamma=3,b_gamma=0.1,
                              a_phi = 4, b_phi = 4,        
                              a_f = 4, b_f = 4,coeff_var = 0.02,
                              cond_order_beta = F, cond_order_gamma = F,cond_order_phi = F,
                              type_warp='lin',WARP = 'Telesca',PRINT = T, HOW_MANY = 10000, 
                              pre_smoot = F, nburn_pre_smoot = 0, niter_pre_smoot = 0,gamma_reg = FALSE,gamma_shrinkage=FALSE)
{
  
  # The function makes the alignment of the curves that takes as input
  
  # There are two types of warping available that are 'Telesca', the one proposed by Telesca 
  # and 'Gamma-Normalization' that is the one proposed by our 
  
  # Two types of smoothing are available regaring the gamma_gi variables
  # The base one is the one with a diagonal covariance matrix (default value)
  # If one sets gamma_reg = TRUE then the prior for gamma_gi becomes the same of the beta_g term
  # gamma_shrinkage is still a working progress for a more complex smoothing component
  
  
  # One can also apply a pre smooting procedure for denoyse the data with the function JARA_smooting_groups
  # Set pre_smoot = TRUE to apply a smoothing of the curves before align them 
  
   
  
  if(WARP != 'Telesca' & WARP != 'Gamma-Normalization'){
    print('WARP specificato non valido')
    return(NA)
  }
  
  if(gamma_reg & gamma_shrinkage){
    print("Only one prior can be assigned to the indvidual parameters")
    return(NA)
  }
  
  n_obs_max  = dim(y)[1] 
  n_patients = dim(y)[2] # y has patients in the columns
  
  
  y_copy = y 
  
  if(pre_smoot){
    y_pre_smoot = JARA_smooting_groups(y, n_groups = n_groups, n_per_group = n_per_group,
                                       nburn=nburn_pre_smoot,niter=niter_pre_smoot,n_knots_m =n_knots_m,n_knots_gamma=n_knots_gamma,
                                       a_eps = a_eps, b_eps = b_eps ,a_lambda = a_lambda, b_lambda = b_lambda,
                                       a_gamma= a_gamma ,b_gamma= b_gamma,
                                       cond_order_beta = cond_order_beta, cond_order_gamma = cond_order_gamma, pre_smoot = T)
    
    y_copy <- matrix(NaN, nrow = n_obs_max, ncol = n_patients)
    col_index <- 1
    for (group in y_pre_smoot) {
      for (patient in group) {
        y_copy[1:length(patient), col_index] <- patient
        col_index <- col_index + 1
      }
    }
    
  }
  
  y_list = list()     # list of lists of groups of lists of individuals for group
  y_val  = list()     # list of lists of numerical values matrices for group 
  n_obs  = list()
  
  n0 = 1 
  nF = n_per_group[1]
  
  for(g in 1:n_groups){
    
    y_temp = y_copy[,n0:nF]
    y_list_temp = lapply(seq_len(ncol(y_temp)), function(i) na.omit(y_temp[, i]))
    n_obs_temp=sapply(y_list_temp,length)
    y_temp = unlist(y_list_temp)
    y_temp = as.numeric(y_temp)
    
    y_list[[g]] = y_list_temp 
    y_val[[g]]  = y_temp 
    n_obs[[g]]  = n_obs_temp
    
    n0 = nF + 1 
    nF = nF + n_per_group[g+1]
  }
  
  rm(n0)
  rm(nF)
  
  if(cond_order_beta){order_beta = 4}else{order_beta = 3}
  if(cond_order_gamma){order_gamma = 4}else{order_gamma = 3}
  if(type_warp=='lin'){degree_phi = 1
  order_phi =2 
  cond_order_phi=T
  }
  
  if(type_warp=='quad'){degree_phi = 2
  if(cond_order_phi){order_phi = 3}else{order_phi = 2}}
  
  if(type_warp=='cub'){ degree_phi = 3 
  if(cond_order_phi){order_phi = 4}else{order_phi = 3}}
  
  
  p = n_knots_m + order_beta      # dim beta 
  k = n_knots_gamma + order_gamma  # dim gamma_i
  q = n_knots_h + order_phi      # dim phi_i
  Omega = omega_P(p)
  P = omega_P(q)
  
  time_group = list()
  
  for(g in 1:n_groups){
    time_group[[g]] = lapply(1:n_per_group[g], function(i)  seq(0,1, length.out = n_obs[[g]][i]))
  }
  
  # knots and Bspline matrix for common shape function m and individual specific term 
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  knots_gamma = seq(0, 1, length.out = n_knots_gamma+2)[-c(1,n_knots_gamma+2)]
  
  Bm_beta_group  = list()
  Bm_gamma_group = list()
  
  
  for(g in 1:n_groups){
    Bm_beta_group[[g]] <- lapply(1:n_per_group[g], function(i) bs(time_group[[g]][[i]], knots = knots_m, intercept = cond_order_beta) )
    Bm_gamma_group[[g]] <- lapply(1:n_per_group[g], function(i) bs(time_group[[g]][[i]], knots = knots_gamma, intercept = cond_order_gamma))
  }
  
  # knots and Bspline matrix for warping functions h_i
  knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
  
  Bh_group = list() 
  
  for(g in 1:n_groups){
    Bh_group[[g]] <- lapply(1:n_per_group[g], function(i) bs(time_group[[g]][[i]], knots = knots_h, intercept = cond_order_phi,degree= degree_phi))
  }
  
  if(WARP == 'Telesca'){
    # Upsilon
    nu <- c(rep(0,order_phi),seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)],rep(1,order_phi))
    Upsilon <- (nu[order_phi] - nu[1])/(order_phi-1)
    for(i in 1:(n_knots_h+order_phi-1)){
      Upsilon[i+1] <- (nu[i+order_phi] - nu[i+1])/(order_phi-1) + Upsilon[i]
    }
  }
  
  
  bigP <- Matrix::bdiag(replicate(n_patients, P, simplify = FALSE))
  
  beta_0 <- rep(0,p)
  gamma_0<- rep(0,k)
  
  tune <- .005
  accepts = list()
  for(g in 1:n_groups){
    accepts[[g]] =  matrix(0, nrow=n_per_group[g], ncol=q-2)
  }
  
  
  #----- Save Structures -----#
  
  nrun <- nburn + niter
  
  #-- Sigma epsilon --#
  sigma_eps_save <- numeric(nrun)
  sigma_eps_save[1] <- 1/rgamma(1,a_eps,rate= b_eps) 
  
  #-- Lambda --#
  lambda_save <- numeric(nrun)
  lambda_save[1] <-1/rgamma(1,a_lambda,rate= b_lambda)  
  
  #-- Beta --#
  beta_save_group <- list()
  for(g in 1:n_groups){
    beta_save_group[[g]] = matrix(NA, nrow = nrun, ncol = p)
    beta_save_group[[g]][1,] <- mvrnorm(n = 1, mu = beta_0, Sigma = solve(Omega/lambda_save[1]) )
  }
  
  X <- lapply(1:n_groups, function(g) stack_Matrix(Bm_beta_group[[g]]))
  
  if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
    #-- Sigma gamma --#
    sigma_gamma_save <- numeric(nrun)
    sigma_gamma_save[1] <- 1/rgamma(1,a_gamma,rate= b_gamma) 
    sigma_gamma_matrix=diag(k)*sigma_gamma_save[1]
    inv_sigma_gamma_matrix = solve(sigma_gamma_matrix)
  }
  if(gamma_reg){
    Omega_gamma = omega_P(k)
    lambda_gamma_save =  numeric(nrun)
    lambda_gamma_save[1] <-1/rgamma(1,a_lambda_gamma,rate= b_lambda_gamma)  
    inv_sigma_gamma_matrix = Omega_gamma/lambda_gamma_save[1]
    sigma_gamma_matrix = solve(inv_sigma_gamma_matrix)
  }
  if(gamma_shrinkage){}
  
  #-- Gamma --#
  gamma_save_group=list() 
  
  for(g in 1:n_groups){
    gamma_save_gruppo_G = list()
    for( i in 1:n_per_group[g]){
      gamma_save_ind = matrix(NA,nrun,k)
      gamma_save_ind[1,] = mvrnorm(n = 1, mu = gamma_0, Sigma = sigma_gamma_matrix )
      gamma_save_gruppo_G[[i]] = gamma_save_ind
    }
    gamma_save_group[[g]] = gamma_save_gruppo_G
  }
  
  rm(gamma_save_gruppo_G)
  rm(gamma_save_ind)
  
  if(WARP == 'Telesca'){
    
    #-- Sigma Phi --#
    sigma_phi_save <- numeric(nrun)
    sigma_phi_save[1] <- 1/rgamma(1,a_phi,rate=b_phi)
    
    
    #-- Phi --#
    phi_save_group = list()   # list of length = nrun where each element is a list per group whose elements
    # are matrices whose rows are patients in that group 
    
    for(g in 1:n_groups){
      phi_save_GROUP_G = list()
      for(i in 1:n_per_group[g]){
        phi_save_ind = matrix(NA,nrun,q)
        phi_save_ind[1,] = Upsilon
        phi_save_GROUP_G[[i]] = phi_save_ind
      }
      phi_save_group[[g]] = phi_save_GROUP_G
    }
    
    phi = list()
    for(g in 1:n_groups){
      phi[[g]] <- lapply(1:n_per_group[g], function(i) Upsilon)
    }
    
    rm(phi_save_GROUP_G)
    rm(phi_save_ind)
    
  }
  if(WARP == 'Gamma-Normalization'){
    
    #-- Csi and phi --#
    csi_save_group = list()
    phi_save_group = list() 
    b_prop_save_group = list()
    
    for(g in 1:n_groups){
      phi_save_GROUP_G = list()
      csi_save_GROUP_G = list()
      b_prop_save_group_G = list()
      for(i in 1:n_per_group[g]){
        phi_save_ind = matrix(NA,nrun,q)
        csi_save_ind = matrix(NA,nrun,q)
        b_prop_save_ind = matrix(NA,nrun-1,q-1)
        
        for(j in 1:q){
          if(j==1){csi_save_ind[1,j] = 0}else{csi_save_ind[1,j] = rgamma(n = 1,a_f[j-1],rate= b_f)}
        }
        
        
        csi_save_GROUP_G[[i]] = csi_save_ind
        phi_save_ind[1,] = cumsum(csi_save_ind[1,])/sum(csi_save_ind[1,])
        phi_save_GROUP_G[[i]] = phi_save_ind
        b_prop_save_group_G[[i]] = b_prop_save_ind
      }
      
      phi_save_group[[g]] = phi_save_GROUP_G
      csi_save_group[[g]] = csi_save_GROUP_G
      b_prop_save_group[[g]] = b_prop_save_group_G
    }
    
    rm(phi_save_GROUP_G)
    rm(csi_save_GROUP_G)
    rm(b_prop_save_group_G)
    rm(phi_save_ind)
    rm(csi_save_ind)
    rm(b_prop_save_ind)
    
    phi = list()
    csi = list()
    for(g in 1:n_groups){
      phi[[g]] <- lapply(1:n_per_group[g], function(i) phi_save_group[[g]][[i]][1,])
      csi[[g]] <- lapply(1:n_per_group[g], function(i) csi_save_group[[g]][[i]][1,])
    }
    
  }
  
  #-- h --#
  h = list()
  
  for(g in 1:n_groups){
    h[[g]] <- lapply(1:n_per_group[g], function(i) time_group[[g]][[i]])
  }
  
  
  #----- MCMC Loop -----#
  
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(iter in 2:nrun){
    
    setTxtProgressBar(bar, iter)
    
    #-- Update Phi --#
    if(WARP == 'Telesca'){
      
      for(g in 1:n_groups){
        
        for (i in 1:n_per_group[g]){
          
          phi_old = phi[[g]][[i]]
          phi_new=phi_old
          
          for (j in 2:(q-1)){
            
            log_densita_old = pi_phi_Tel(phi_old,beta_save_group[[g]][iter-1,], gamma_save_group[[g]][[i]][iter-1,], 
                                         sigma_eps_save[iter-1], P,sigma_phi_save[iter-1], Upsilon,Bh_group[[g]][[i]],y_list[[g]][[i]],cond_order_beta,cond_order_gamma)
            
            phi_new[j] <- Arma_runif(max(phi_new[j] - tune, phi_new[j-1]),min(phi_new[j] + tune, phi_new[j+1]))
            
            log_densita_new=pi_phi_Tel(phi_new,beta_save_group[[g]][iter-1,], gamma_save_group[[g]][[i]][iter-1,], 
                                       sigma_eps_save[iter-1], P,sigma_phi_save[iter-1], Upsilon,Bh_group[[g]][[i]],y_list[[g]][[i]],cond_order_beta,cond_order_gamma)
            
            alpha = min(0, (log_densita_new-log_densita_old))
            
            u = Arma_runif(0, 1)
            
            if(u<exp(alpha)){#accept 
              accepts[[g]][i, j-1]= accepts[[g]][i, j-1]+1} 
            else{ phi_new[j]= phi_old[j]}
          }
          phi[[g]][[i]] = phi_new
          phi_save_group[[g]][[i]][iter,] = phi_new
        }
      }
      
      #-- Update sigma_phi --#
      a_star = a_phi + 0.5*q*n_patients
      diff_term= lapply(1:n_groups, function(g) lapply(1:n_per_group[g],function(i){phi_save_group[[g]][[i]][iter,]-Upsilon}))  
      prod_term=  lapply(1:n_groups, function(g) lapply(diff_term[[g]],function(v){t(v)%*%P%*%v})) 
      somma=sum(unlist(prod_term))
      b_star=b_phi+0.5*somma
      sigma_phi_save[iter]=1/rgamma(1,a_star,rate=b_star)
      
    }
    
    if(WARP=='Gamma-Normalization'){
      
      for(g in 1:n_groups){
        
        for (i in 1:n_per_group[g]){
          
          phi_old = phi[[g]][[i]]
          phi_new = phi_old
          
          csi_old = csi[[g]][[i]]
          csi_new = csi_old
          
          for (j in 2:(q-1)){
            
            log_densita_old = pi_phi_GN(phi_old,beta_save_group[[g]][iter-1,], gamma_save_group[[g]][[i]][iter-1,], 
                                        sigma_eps_save[iter-1],csi_old[j],a_f[j-1],b_f,Bh_group[[g]][[i]],y_list[[g]][[i]],cond_order_beta,cond_order_gamma)
            
            a_prop= 1/coeff_var^2
            b_prop= a_prop/csi_old[j]
            b_prop_save_group[[g]][[i]][iter-1,j] = b_prop
            
            csi_new[j]=rgamma(n=1,a_prop,b_prop) 
            phi_new<-cumsum(csi_new)/sum(csi_new)
            
            log_densita_new = pi_phi_GN(phi_new,beta_save_group[[g]][iter-1,], gamma_save_group[[g]][[i]][iter-1,], 
                                        sigma_eps_save[iter-1],csi_new[j],a_f[j-1],b_f,Bh_group[[g]][[i]],y_list[[g]][[i]],cond_order_beta,cond_order_gamma)
            
            quoziente=(log_densita_new-log_densita_old)+
              (dgamma(csi_old[j],a_prop,b_prop,log=T)-dgamma(csi_new[j],a_prop,b_prop,log=T))
            
            alpha = min(0, quoziente)
            
            u = Arma_runif(0, 1)
            
            if(PRINT){
              if(iter%%HOW_MANY == 0){
                cat("Patient =", i, "\n")
                cat("Log_dens_new =", log_densita_new, "\n")
                cat("Log_dens_old =", log_densita_old, "\n")
                cat("alpha =", alpha, "\n")
                cat("Accepted =", u <= exp(alpha), "\n")
                cat("\n")
              }
            }
            
            if(u<exp(alpha)){#accept 
              accepts[[g]][i, j-1]= accepts[[g]][i, j-1]+1} 
            else{ 
              phi_new[j] = phi_old[j]
              csi_new[j] = csi_old[j]}
          }
          phi[[g]][[i]] = phi_new
          phi_save_group[[g]][[i]][iter,] = phi_new
          csi[[g]][[i]] = csi_new
          csi_save_group[[g]][[i]][iter,] = csi_new
        }
      }
    }
    
    #-- Update h  --#
    
    h <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) as.numeric(Bh_group[[g]][[i]] %*% phi_save_group[[g]][[i]][iter,])) )
    
    #-- Update Bm_beta and Bm_gamma --# 
    
    Bm_beta_group  = lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) bs(h[[g]][[i]], knots = knots_m, intercept = cond_order_beta) ) )
    Bm_gamma_group = lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) bs(h[[g]][[i]], knots = knots_gamma, intercept = cond_order_gamma)) )
    
    
    #-- Update X --#
    
    X <-  lapply(1:n_groups, function(g) stack_Matrix(Bm_beta_group[[g]]))
    
    #-- Update Beta --#
    
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = lapply(1:n_groups, function(g) inv_sigma_beta + crossprod(X[[g]])/(sigma_eps_save[iter-1]) )    
    
    V_beta = lapply(1:n_groups, function(g)   armaInv(inv_V_beta[[g]])  )
    
    vec_list= lapply(1:n_groups, function(g)  lapply(1:n_per_group[g], function(i) gamma_save_group[[g]][[i]][iter-1,]) )    
    G = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { matrice %*% vettore}, Bm_gamma_group[[g]], vec_list[[g]]))   
    C_gamma<- lapply(1:n_groups, function(g) rep(do.call(rbind, G[[g]]),1)) 
    
    m_beta = lapply(1:n_groups, function(g)   V_beta[[g]]%*%(t(X[[g]]) %*% ((y_val[[g]]-C_gamma[[g]])/sigma_eps_save[iter-1]) + inv_sigma_beta %*% beta_0 )  ) 
    Beta_udpate =  lapply(1:n_groups, function(g) Arma_mvrnorm(m_beta[[g]], V_beta[[g]]) ) 
    
    aggiorna_matrice <- function(matrice, vettore) {
      matrice[iter, ] <- vettore
      return(matrice)
    }
    
    beta_save_group <- mapply(aggiorna_matrice,  beta_save_group,  Beta_udpate , SIMPLIFY = FALSE)
    
    #-- Update lambda --#
    
    a_star = 0.5*p*n_groups + a_lambda
    beta_somma = lapply(1:n_groups, function(g) t(beta_save_group[[g]][iter,])%*% Omega %*%beta_save_group[[g]][iter,] )
    sum_term = as.numeric( sum(unlist(beta_somma)) )
    b_star = 0.5 * sum_term +  b_lambda
    lambda_save[iter]<-1/rgamma(n=1, shape=a_star, rate=b_star)
    
    
    #-- Update Gamma --#
    if(gamma_shrinkage==FALSE){
      # The update is the same for the case with or witjot reguralization, onlt inv-sigma matrix changes but it is updated separately
      prod_X_tilde =  lapply(1:n_groups, function(g) lapply(Bm_gamma_group[[g]],function(matrice){crossprod(matrice)/sigma_eps_save[iter-1]})   )   
      inv_V_gamma  =  lapply(1:n_groups, function(g) lapply(prod_X_tilde[[g]],function(matrice){matrice + inv_sigma_gamma_matrix})  ) 
      V_gamma      =  lapply(1:n_groups, function(g) lapply(inv_V_gamma[[g]], function(matrice){armaInv(matrice)}) ) 
      
      c_tilde   = lapply(1:n_groups, function(g) lapply(Bm_beta_group[[g]],function(matrice){matrice%*%beta_save_group[[g]][iter,]})  ) 
      diff_term = lapply(1:n_groups, function(g) Map('-',y_list[[g]],c_tilde[[g]]) )
      prod_term = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { t(matrice) %*% vettore}, Bm_gamma_group[[g]], diff_term[[g]]) )
      m_gamma   = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { (1/sigma_eps_save[iter-1])*matrice %*% vettore}, V_gamma[[g]], prod_term[[g]]) ) 
      
      G = lapply(1:n_groups, function(g) Map(function(matrice, vettore) {mvrnorm(n = 1, mu = vettore, Sigma = matrice )},
                                             V_gamma[[g]], m_gamma[[g]]) ) 
      
      gamma_save_group <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g],
                                                                function(i) aggiorna_matrice(gamma_save_group[[g]][[i]],  G[[g]][[i]] ) ) ) 
    }
    
    # Update gamma variance parameters 
    
    if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
      #-- Update Sigma gamma --#
      
      a_star=0.5*n_patients*k+a_gamma
      
      vettore_somma = lapply(1:n_groups, function(g)  lapply(1:n_per_group[g],function(i){crossprod(gamma_save_group[[g]][[i]][iter,])}) )
      b_star = b_gamma + 0.5 * sum(unlist(vettore_somma))
      
      sigma_gamma_save[iter]=1/rgamma(1,a_star,rate=b_star)
      
      inv_sigma_gamma_matrix=diag(k)/sigma_gamma_save[iter]
    }
    if(gamma_reg){
      #-- Update lambda_gamma --#  
      a_star=0.5*n_patients*k+a_gamma
      vettore_somma = lapply(1:n_groups, function(g)  lapply(1:n_per_group[g],function(i){ t(gamma_save_group[[g]][[i]][iter,])%*% Omega_gamma %*% gamma_save_group[[g]][[i]][iter,]}) )
      b_star = b_gamma + 0.5 * sum(unlist(vettore_somma))
      lambda_gamma_save[iter] = 1/rgamma(n=1, shape=a_star, rate=b_star)
      inv_sigma_gamma_matrix = Omega_gamma/lambda_gamma_save[iter]
    }
    
    
    #-- Update sigma eps --#
    
    a_star = a_eps + 0.5*sum(unlist(n_obs))
    
    term_beta = lapply(1:n_groups, function(g) lapply(Bm_beta_group[[g]],function(bm){bm%*%beta_save_group[[g]][iter,]}) )
    
    #vec_list è gia calcolato da prima 
    #vec_list= lapply(1:n_groups, function(g)  lapply(1:n_per_group[g], function(i) gamma_save_group[[g]][[i]][iter-1,]) ) 
    term_gamma = lapply(1:n_groups, function(g) Map(function(matrice, vettore) { matrice %*% vettore}, Bm_gamma_group[[g]], vec_list[[g]]) )
    
    m_i = lapply(1:n_groups, function(g) Map("+", term_beta[[g]], term_gamma[[g]]) )
    diff_term = lapply(1:n_groups, function(g) Map("-", y_list[[g]], m_i[[g]]) ) 
    prod_term = lapply(1:n_groups, function(g) lapply(diff_term[[g]],function(v){crossprod(v)}) )
    somma=0.5*sum(unlist(prod_term))
    b_star=somma+b_eps
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    if(iter%%1000 == 0){ print(iter)}
    
    
  }
  close(bar)
  
  #-- Acceptance rates --# 
  accepts = lapply(1:n_groups, function(g) accepts[[g]] / nrun )
  
  #-- Posterior paramameters estimation --#
  
  beta_post_group  = lapply( 1: n_groups, function(g) colMeans(beta_save_group[[g]][-c(1:nburn),]) )
  
  gamma_post_group = lapply(1: n_groups, function(g) lapply(1:n_per_group[g], 
                                                            function(i) colMeans(gamma_save_group[[g]][[i]][-c(1:nburn),]) ))
  
  lambda_post <- mean(lambda_save[-c(1:nburn)])
  sigma_eps_post<-mean(sigma_eps_save[-c(1:nburn)])
  if(gamma_reg==FALSE & gamma_shrinkage==FALSE){sigma_gamma_post<-mean(sigma_gamma_save[-c(1:nburn)])}
  if(gamma_reg){lambda_gamma_post = mean(lambda_gamma_save[-c(1:nburn)])}
  if(gamma_shrinkage){}
  phi_post_group   = lapply(1: n_groups, function(g) lapply(1:n_per_group[g], 
                                                            function(i) colMeans(phi_save_group[[g]][[i]][-c(1:nburn),]) ))
  
  if(WARP == 'Telesca'){
    sigma_phi_post<-mean(sigma_phi_save[-c(1:nburn)])
  }
  if(WARP == 'Gamma-Normalization'){
    csi_post_group   = lapply(1: n_groups, function(g) lapply(1:n_per_group[g], 
                                                              function(i) colMeans(csi_save_group[[g]][[i]][-c(1:nburn),]) ))
  }
  
  #h_p <- lapply(h_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  h_p = lapply(1:n_groups, function(g) lapply(1:n_per_group[g],function(i) Bh_group[[g]][[i]]%*%phi_post_group[[g]][[i]]) ) 
  
  Bm_post_beta_group  = lapply(1:n_groups, function(g) lapply(h_p[[g]], function(h) bs(h, knots = knots_m, intercept = cond_order_beta)) ) 
  Bm_post_gamma_group = lapply(1:n_groups, function(g) lapply(h_p[[g]], function(h) bs(h, knots = knots_gamma, intercept = cond_order_gamma)) ) 
  
  y_p <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) 
    as.numeric(Bm_post_beta_group[[g]][[i]] %*% beta_post_group[[g]] + Bm_post_gamma_group[[g]][[i]] %*% gamma_post_group[[g]][[i]])) ) 
  y_star <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) interp_spline(h_p[[g]][[i]], y_list[[g]][[i]])) ) 
  y_star_smoot <- lapply(1:n_groups, function(g) lapply(1:n_per_group[g], function(i) interp_spline(h_p[[g]][[i]], y_p[[g]][[i]])) ) 
  
  print('DONE')
  
  if(pre_smoot){
    if(WARP == 'Telesca'){
      if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, sigma_gamma=sigma_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,sigma_phi = sigma_phi_post,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, sigma_gamma_s=sigma_gamma_save,phi_group_s=phi_save_group,sigma_phi_s=sigma_phi_save), n_obs = n_obs,
                      y_star=y_star,y_star_smoot=y_star_smoot,Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group),
                      pre_smoot = list(post = list(y_p = y_pre_smoot))))
      }
      if(gamma_reg){ 
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, lambda_gamma=lambda_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,sigma_phi = sigma_phi_post,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, lambda_gamma_s=lambda_gamma_save,phi_group_s=phi_save_group,sigma_phi_s=sigma_phi_save), n_obs = n_obs,
                      y_star=y_star,y_star_smoot=y_star_smoot,Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group),
                      pre_smoot = list(post = list(y_p = y_pre_smoot))))
      }
      if(gamma_shrinkage){}
    }
    if(WARP == 'Gamma-Normalization'){
      if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, sigma_gamma=sigma_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,csi_group = csi_post_group,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, sigma_gamma_s=sigma_gamma_save,phi_group_s=phi_save_group,csi_group_s=csi_save_group,b_prop_group_s = b_prop_save_group),
                      n_obs = n_obs,y_star=y_star,y_star_smoot=y_star_smoot,
                      Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group),
                      pre_smoot = list(post = list(y_p = y_pre_smoot)) ))
      }
      if(gamma_reg){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, lambda_gamma=lambda_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,csi_group = csi_post_group,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, lambda_gamma_s=lambda_gamma_save,phi_group_s=phi_save_group,csi_group_s=csi_save_group,b_prop_group_s = b_prop_save_group),
                      n_obs = n_obs,y_star=y_star,y_star_smoot=y_star_smoot,
                      Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group),
                      pre_smoot = list(post = list(y_p = y_pre_smoot)) ))
      }
      if(gamma_shrinkage){}
    }
  }
  else{
    if(WARP == 'Telesca'){
      if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, sigma_gamma=sigma_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,sigma_phi = sigma_phi_post,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, sigma_gamma_s=sigma_gamma_save,phi_group_s=phi_save_group,sigma_phi_s=sigma_phi_save), n_obs = n_obs,
                      y_star=y_star,y_star_smoot=y_star_smoot,Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group)))
      }
      if(gamma_reg){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, lambda_gamma=lambda_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,sigma_phi = sigma_phi_post,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, lambda_gamma_s=lambda_gamma_save,phi_group_s=phi_save_group,sigma_phi_s=sigma_phi_save), n_obs = n_obs,
                      y_star=y_star,y_star_smoot=y_star_smoot,Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group)))
      }
      if(gamma_shrinkage){}
    }
    if(WARP == 'Gamma-Normalization'){
      if(gamma_reg==FALSE & gamma_shrinkage==FALSE){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, sigma_gamma=sigma_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,csi_group = csi_post_group,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, sigma_gamma_s=sigma_gamma_save,phi_group_s=phi_save_group,csi_group_s=csi_save_group,b_prop_group_s = b_prop_save_group),
                      n_obs = n_obs,y_star=y_star,y_star_smoot=y_star_smoot,
                      Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group)))
      }
      if(gamma_reg){
        return (list (post=list(beta_group = beta_post_group, lambda = lambda_post,  
                                sigma_eps = sigma_eps_post, lambda_gamma=lambda_gamma_post,
                                gamma_group=gamma_post_group, y=y_p, accepts=accepts, phi_group =phi_post_group ,csi_group = csi_post_group,h=h_p), 
                      full=list(beta_group_s=beta_save_group, lambda_s=lambda_save, sigma_eps_s=sigma_eps_save,
                                gamma_group_s=gamma_save_group, lambda_gamma_s=lambda_gamma_save,phi_group_s=phi_save_group,csi_group_s=csi_save_group,b_prop_group_s = b_prop_save_group),
                      n_obs = n_obs,y_star=y_star,y_star_smoot=y_star_smoot,
                      Bm = list(Bm_beta_group = Bm_post_beta_group, Bm_gamma_group = Bm_post_gamma_group, Bh_group= Bh_group)))
      }
      if(gamma_shrinkage){}
    }
  }
  
}


pre_smoot = F
niter_pre_smoot = 0
nburn_pre_smoot = 0

media_inv_gamma = function(a,b){b/(a-1)}
var_inv_gamma   = function(a){media_inv_gamma(a,b)^2/(a-2)}
media_var_phi_ij    = function(a,j){
  if(j == length(a)){
    media = 1 
    varianza = 0
  }
  else{
    parz  = cumsum(a)[j]
    da_j  = sum(a[-(1:j)])
    total = sum(a)
    media = parz/ total
    varianza = parz * da_j /((total-1)*(total)^2)
  }
  return(list(media=media,varianza=varianza))
}
