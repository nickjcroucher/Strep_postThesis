# save some functions

build_K <- function(beta, m) {
  K <- matrix(NA, nrow = 2, ncol = 2)
  K[1, 1] <- beta*m[1, 1]*T1+beta*m[1, 2]*G
  K[1, 2] <- beta*m[1, 2]*T2
  K[2, 1] <- beta*m[2, 1]*T1+beta*m[2, 2]*G
  K[2, 2] <- beta*m[2, 2]*T2
  return(K)
}

# edit R0 as optional parameter for elastiticy analysis for sigmas & deltas
compute_R0 <- function(pars,
                       
                       sigma1_group1 = NULL,
                       sigma1_group2 = NULL,
                       
                       delta_group1 = NULL,
                       delta_group2 = NULL
){
  
  # sigma & deltas are editable for each age groups!
  s1 <- (if(!is.null(sigma1_group1)) sigma1_group1 else
    pars$sigma_1)
  s2 <- (if(!is.null(sigma1_group2)) sigma1_group2 else
    pars$sigma_1)
  
  d1 <- (if(!is.null(delta_group1)) delta_group1 else
    pars$delta_1)
  d2 <- (if(!is.null(delta_group2)) delta_group2 else
    pars$delta_2)
  
  rA1 <- d1+s1+pars$mu_0_1+pars$alpha
  rA2 <- d2+s2+pars$mu_0_2
  rD1 <- pars$sigma_2+pars$mu_1+pars$mu_0_1+pars$alpha
  rD2 <- pars$sigma_2+pars$mu_1+pars$mu_0_2
  
  T1 <- (rD1+d1)/(rA1*rD1)
  T2 <- (rD2+d2)/(rA2*rD2)
  G <- (pars$alpha/rA1)*(T2+d1/(rD1*rD2))
  
  K <- matrix(NA, nrow = 2, ncol = 2)
  K[1, 1] <- pars$beta_0*m[1, 1]*T1+pars$beta_0*m[1, 2]*G
  K[1, 2] <- pars$beta_0*m[1, 2]*T2
  K[2, 1] <- pars$beta_0*m[2, 1]*T1+pars$beta_0*m[2, 2]*G
  K[2, 2] <- pars$beta_0*m[2, 2]*T2
  
  a <- K[1, 1]; b <- K[1, 2]
  c <- K[2, 1]; d <- K[2, 2]
  
  # algebraic R0
  return(((a+d)+sqrt((a-d)^2+4*b*c))/2)
}

elasticity_param <- function(param_name,
                             pars,
                             delta = 0.001){
  R0_base <- compute_R0(pars)
  pars_high <- pars
  pars_high[[param_name]] <- pars[[param_name]]*(1 + delta)
  R0_high <- compute_R0(pars_high)
  
  # Chitnis et al. (2008)
  return(((R0_high-R0_base)/(pars[[param_name]]*delta))*
           (pars[[param_name]]/R0_base))
}

elasticity_sigma <- function(group,
                             pars,
                             delta = 0.001){
  sigma_base <- pars$sigma_1
  sigma_high <- sigma_base*(1 + delta)
  
  R0_base <- compute_R0(pars,
                        sigma1_group1 = sigma_base,
                        sigma1_group2 = sigma_base)
  
  # perturb sigma1 just in children OR adults
  R0_high <- if (group == "group1"){
    compute_R0(pars,
               sigma1_group1 = sigma_high,
               sigma1_group2 = sigma_base)
  } else {
    compute_R0(pars,
               sigma1_group1 = sigma_base,
               sigma1_group2 = sigma_high)
  }
  
  # Chitnis et al. (2008)
  return(((R0_high-R0_base)/(sigma_base*delta))*
           (sigma_base/R0_base))
}
