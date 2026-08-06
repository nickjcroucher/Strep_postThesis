# See https://mrc-ide.github.io/mcstate/articles/nested_sir_models.html

burnin_days <- 0
# 
# ll_nbinom <- function(data, model, kappa, exp_noise) {
#   # if (is.na(data)) {
#   #   return(numeric(length(model)))
#   # }
# 
#   data_clean <- ifelse(is.na(data), 0, data)
#   mu <- model + rexp(length(model), rate = exp_noise)
#   dnbinom(data_clean, kappa, mu = mu, log = TRUE)
# }
# 
# case_compare <- function(state, observed, pars = NULL) {
#   exp_noise <- 1e6
#   n <- ncol(state)
#   # kappa_1 <- 5
# 
#   # sir_model$info()$index$n_AD_weekly
#   model_1 <- state[7, , drop = TRUE]
#   model_2 <- state[8, , drop = TRUE]
# 
#   if (is.na(observed$count_s1_1)) {
#     ll_1 <- ll_nbinom(data = 0,
#                          model = model_1,
#                          kappa = pars$kappa_1,
#                          exp_noise = exp_noise)
#   } else {
#     ll_1 <- ll_nbinom(data = observed$count_s1_1,
#                          model = model_1,
#                          kappa = pars$kappa_1,
#                          exp_noise = exp_noise)
#   }
# 
#   if (is.na(observed$count_s1_2)) {
#     ll_2 <- ll_nbinom(data = 0,
#                          model = model_2,
#                          kappa = pars$kappa_1,
#                          exp_noise = exp_noise)
#   } else {
#     ll_2 <- ll_nbinom(data = observed$count_s1_2,
#                          model = model_2,
#                          kappa = pars$kappa_1,
#                          exp_noise = exp_noise)
#   }
# 
#   ll <- ll_1 + ll_2
#   return(ll)
# }

case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6

  # sir_model$info()$index$n_AD_weekly
  model_1 <- state[7, , drop = TRUE]
  model_2 <- state[8, , drop = TRUE]

  if (is.na(observed$count_s1_1)) {
    ll_1 <- dpois(x = 0,
                  lambda = model_1 + rexp(ncol(state), exp_noise),
                  log = TRUE
    )
  } else {
    ll_1 <- dpois(x = observed$count_s1_1,
                  lambda = model_1 + rexp(ncol(state), exp_noise),
                  log = TRUE
    )
  }

  if (is.na(observed$count_s1_2)) {
    ll_2 <- dpois(x = 0,
                  lambda = model_2 + rexp(ncol(state), exp_noise),
                  log = TRUE
    )
  } else {
    ll_2 <- dpois(x = observed$count_s1_2,
                  lambda = model_2 + rexp(ncol(state), exp_noise),
                  log = TRUE
    )
  }

  ll <- ll_1 + ll_2
  return(ll)
}


# generate index function
index_fun <- function(info){
  if (is.null(info$index)){
    info <- info[[1]]
  }
  list(run = unlist(info$index),
       state = unlist(info$index))
}

# That transform function
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/tests/testthat/test-if2.R#L40
# https://github.com/mrc-ide/mcstate/issues/184
parameter_transform <- function(t_norm) {
  library(socialmixr)
  age.limits = c(0, 15)
  N_age <- length(age.limits)
  
  contact_2_demographic <- suppressMessages(
    socialmixr::contact_matrix(polymod,
                               countries = "United Kingdom",
                               age.limits = age.limits,
                               symmetric = TRUE
    ))
  
  transmission <- contact_2_demographic$matrix /
    rep(contact_2_demographic$demography$population,
        each = ncol(contact_2_demographic$matrix))
  t_norm <- transmission/max(transmission)
  
  transform <- function(pars) {
    # re-define pars with pars that I really wanna fit only
    log_A_ini <- pars[["log_A_ini"]] # pars[paste0("log_A_ini", 1:2)]
    phi <- pars[["phi"]]
    
    time_shift_1 <- pars[["time_shift_1"]]
    beta_0 <- pars[["beta_0"]]
    beta_1 <- pars[["beta_1"]]
    # beta_diff <- pars[["beta_diff"]]
    vacc <- pars[["vacc"]]
    
    
    log_delta1 <- pars[["log_delta1"]]
    # rho <- pars[["rho"]]
    log_delta2 <- pars[["log_delta2"]]
    # sigma_1 <- pars[["sigma_1"]]
    omega <- pars[["omega"]]
    # kappa_1 <- pars[["kappa_1"]]
    
    pars <- list(log_A_ini = log_A_ini,
                 phi = phi,
                 time_shift_1 = time_shift_1,
                 beta_0 = beta_0,
                 beta_1 = beta_1,
                 # beta_diff = beta_diff,
                 vacc = vacc,
                 log_delta1 = log_delta1,
                 # rho = rho
                 log_delta2 = log_delta2,
                 # sigma_1 = sigma_1,
                 omega = omega
                 # kappa_1 = kappa_1
    )
    
    pars$N_ini <-  contact_2_demographic$demography$population
    pars$m <- t_norm
    
    pars
  }
  
  transform
}
# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/transformation
transform <- parameter_transform(t_norm)

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(
      mcstate::pmcmc_parameter("log_A_ini", 0.55, min = 0, max = 1,
                               prior = priors$log_A_ini),
      mcstate::pmcmc_parameter("phi", 0.94, min = (0), max = 2,
                               prior = priors$phi),
      mcstate::pmcmc_parameter("time_shift_1", 0.2, min = (0), max = 0.5, # previously (-10, 1)
                               prior = priors$time_shifts),
      mcstate::pmcmc_parameter("beta_0", 0.05959, min = 0, max = 0.5, # max based on 1/values; worst case increased to 5x
                               prior = priors$beta_0),
      mcstate::pmcmc_parameter("beta_1", 0.2, min = 0, max = 1,
                               prior = priors$betas),
      # mcstate::pmcmc_parameter("beta_diff", 0.8, min = 0, max = 1,
      #                          prior = priors$beta_diff),
      mcstate::pmcmc_parameter("vacc", 0.0001, min = 0, max = 1,
                               prior = priors$vacc),
      mcstate::pmcmc_parameter("log_delta1", (-4), min = (-8), max = -2, #(-3.8), min = (-5), max = -2, #-0.03196764, # log10(1/UK_calibration_kids) for delta1 = 1
                               prior = priors$log_delta1),
      mcstate::pmcmc_parameter("log_delta2", (-3.5), min = (-8), max = -2, #(-3.8), min = (-5), max = -2, #-0.03196764, # log10(1/UK_calibration_kids) for delta1 = 1
                               prior = priors$log_delta2),
      # mcstate::pmcmc_parameter("sigma_1", 0.063, min = 0, max = 1,
      #                          prior = priors$sigmas),
      mcstate::pmcmc_parameter("omega", 2e-4, min = 0, max = 1,
                               prior = priors$omega)
    #   mcstate::pmcmc_parameter("kappa_1", 10, min = 0,
    #                            prior = priors$kappas)
    ),
    proposal = proposal,
    transform = transform
  )
}

prepare_priors <- function(pars) {
  priors <- list()
  
  priors$log_A_ini <- function(s) {
    dnorm(s, mean = 0.54, sd = 0.5, log = TRUE)
    # dnorm(s, mean = 0.54, sd = 0.05, log = TRUE)
  }
  priors$phi <- function(s) {
    dnorm(s, mean = 0.94, sd = 0.1, log = TRUE) # previously 0.15, 0.05
    # stabledist::dstable(s, alpha = 2, beta = 0, gamma = 0.3, delta = 6, log = TRUE) # previously 0.5
  }
  priors$time_shifts <- function(s) {
    dgamma(s, shape=2, scale=0.08, log=TRUE)
    # stabledist::dstable(s, alpha = 2, beta = 0, gamma = 0.5, delta = -5, log = TRUE)
  }
  priors$beta_0 <- function(s) {
    dbeta(s, 3, 100, log = TRUE)
    # dbeta(s, 30, 500, log = TRUE) # or dbeta(s, 20, 400, log = TRUE)
  }
  priors$betas <- function(s) {
    dbeta(s, 3.5, 10, log = TRUE) # more relaxed dbeta(s, 2, 15
    # dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$beta_diff <- function(s) {
    dnorm(s, mean = 0.8, sd = 0.1, log = TRUE) # more relaxed dbeta(s, 2, 15
    # dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$vacc <- function(s) {
    dbeta(s, 1, 1.5, log = TRUE)
    # dnorm(s, mean = 1e-4, sd = 3e-5, log = TRUE)
  }
  priors$log_delta1 <- function(s) {
    dnorm(s, mean = -4, sd = 0.9, log = TRUE)
    # dnorm(s, mean = -4, sd = 0.12, log = TRUE)
  }
  priors$log_delta2 <- function(s) {
    dnorm(s, mean = -4, sd = 0.9, log = TRUE)
    # dnorm(s, mean = -3.42, sd = 0.1, log = TRUE)  # mean=-3.5, sd=0.15
  }
  priors$rho <- function(s) {
    # dnorm(s, mean = 0.8, sd = 0.5, log = TRUE)
    # dgamma(s, shape=4, scale=0.15, log=TRUE) # avoid Cauchy
    stabledist::dstable(s, alpha = 2, beta = 0, gamma = 0.5, delta = -4.5, log = TRUE) # alpha = 1 = Cauchy (?)
    # dunif(s, min = -1, max = 2, log = TRUE)
  }
  priors$sigmas <- function(s) {
    dnorm(s, mean = 0.063, sd = 0.003, log = TRUE)
  }
  priors$omega <- function(s) {
    dunif(s, min = 0, max = 0.5, log = TRUE)
    # dgamma(s, shape = 4, scale = 5e-5, log = TRUE)
  }
  priors$kappas <- function(s) {
    dgamma(s, shape=7,scale=1, log = TRUE)
    # stabledist::dstable(s, alpha = 2, beta = 0, gamma = 1.5, delta = 5, log = TRUE)
    # stabledist::dstable(s, alpha = 2, beta = 0, gamma = 3, delta = 10, log = TRUE)
    # dunif(s, min = 0, max = 100, log = TRUE)
  }
  
  priors
}


pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = round(n_steps*0.5), thin = NULL)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  mcmc1 <- coda::as.mcmc(cbind(pmcmc_result$probabilities, pmcmc_result$pars))
  mcmc1
}

ess_calculation <- function(mcmc1){
  # compile par names & generate switch
  par_names <- colnames(mcmc1)
  
  ess_values <- sapply(par_names, function(p){
    trace <- mcmc1[, p]
    if (var(trace) == 0) {
      warning(sprintf("Parameter '%s' has zero variance. ESS set to NA.", p))
      return(NA)
    } else {
      return(coda::effectiveSize(trace))
    }
  })
  acceptance_rate = 1 - coda::rejectionRate(mcmc1)
  
  list(
    ess = ess_values,
    acceptance_rate = acceptance_rate
  )
}

pmcmc_trace <- function(mcmc1) {
  plot(mcmc1) # to save the figures into pdf
}

################################################################################
# Tuning functions
tuning_pmcmc_further_process <- function(n_steps, tune_pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(tune_pmcmc_result,
                                          burnin = round(n_steps*0.5),
                                          thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  tune_pmcmc_result <- coda::as.mcmc(cbind(processed_chains$probabilities,
                                           processed_chains$pars))
  tune_pmcmc_result
}

################################################################################
# MCMC Diagnostics
# 1. Gelman-Rubin Diagnostic
# https://cran.r-project.org/web/packages/coda/coda.pdf

diag_init_gelman_rubin <- function(tune_pmcmc_result){
  n_chains <- 4 # tune_control$n_chains
  n_samples <- nrow(tune_pmcmc_result$pars)/n_chains
  
  # Split the parameter samples and probabilities by chains
  chains <- lapply(1:n_chains, function(i) {
    start <- (i - 1) * n_samples + 1
    end <- i * n_samples
    list(
      pars = tune_pmcmc_result$pars[start:end, ],
      probabilities = tune_pmcmc_result$probabilities[start:end, ]
    )
  })
  
  # Convert chains to mcmc objects
  mcmc_chains <- lapply(chains, function(chain) {
    coda::as.mcmc(cbind(chain$probabilities, chain$pars))
  })
  
  # Combine chains into a list
  mcmc_chains_list <- do.call(coda::mcmc.list, mcmc_chains)
  mcmc_chains_list
}

diag_cov_mtx <- function(mcmc_chains_list) {
  # print("Covariance matrix of mcmc2")
  cov(as.matrix(mcmc_chains_list))
}

diag_gelman_rubin <- function(mcmc_chains_list) {
  # print("Gelman-Rubin diagnostic")
  gelman_plot <- coda::gelman.plot(mcmc_chains_list,
                                   bin.width = 10,
                                   max.bins = 50,
                                   confidence = 0.95,
                                   transform = FALSE,
                                   autoburnin=TRUE,
                                   auto.layout = TRUE)
  
  coda::gelman.diag(mcmc_chains_list,
                    confidence = 0.95,
                    transform=FALSE,
                    autoburnin=TRUE,
                    multivariate=F)
  
}

# 2. Autocorrelation plots
diag_aucorr <- function(mcmc2){
  for (name in colnames(mcmc2)){
    print(coda::acfplot(mcmc2[, name], main = name))
  }
}

################################################################################
# Particle samples (adapted from Lilith's)
observe_pois <- function(lambda) {
  n_par <- nrow(lambda)
  n_obs <- ncol(lambda)
  ret <- vapply(seq_len(n_par), function(i) {
    rpois(n_obs, lambda[i, ])}, numeric(n_obs))
  t(ret)
}

observe <- function(pmcmc_samples) {
  
  state <- pmcmc_samples$trajectories$state
  pars <- apply(pmcmc_samples$pars, MARGIN = 1, pmcmc_samples$predict$transform)
  time <- pmcmc_samples$trajectories$time
  
  ## extract model outputs
  model_all <- (state[7, , , drop = TRUE]+state[8, , , drop = TRUE])
  model_child <- state[7, , , drop = TRUE]
  
  observed <- list()
  observed$cases_child <- observe_pois(model_child)
  
  abind::abind(c(list(state), observed), along = 1)
}

plot_states <- function(state, data) {
  col <- grey(0.3, 0.1)
  # model state refers to n_AD_weekly (not the D compartment)
  matplot(data$yearWeek, t((state["n_AD1_weekly", , -1]+state["n_AD2_weekly", , -1])),
          type = "l", lty = 1, col = col, ylim = c(0, 41),
          xlab = "", ylab = "Serotype 1 cases")
  # points(data$yearWeek, data$count_serotype, col = 3, pch = 20)
  points(data$yearWeek, (data$count_s1_1+data$count_s1_2), col = 4, type = "l")
  
  matplot(data$yearWeek, xlab = "", t(state["S_tot", , -1]),
          type = "l", lty = 1, col = 2, ylab = "%", ylim = c(0, 6.7e7), yaxt = "n")
  axis(side = 2, at = seq(0, 6e7, length.out = 5),
       labels = seq(0, 100, length.out = 5))
  
  matlines(data$yearWeek, t(state["A_tot", , -1]), lty = 1, col = 1)
  matlines(data$yearWeek, t(state["D_tot", , -1]), lty = 1, col = 4)
  matlines(data$yearWeek, t(state["R_tot", , -1]), lty = 1, col = 5)
  legend("right", bty = "n", fill = 2:4, legend = c("S_tot", "A_tot", "D_tot", "R_tot"))
  # 
  # matplot(data$yearWeek, xlab = "", t(state["I_tot", , -1]),
  #         type = "l", lty = 1, col = 3, ylab = "carriers")
}
