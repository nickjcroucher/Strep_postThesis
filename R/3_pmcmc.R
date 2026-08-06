# 2. Data Fitting ##############################################################
library(mcstate)
library(coda)
library(odin.dust)
library(dust)
library(GGally)
library(socialmixr)

source("global/all_function_allAge.R")
sir_data <- readRDS("raw_data/pmcmc_data_week_allAge_ser1_test_2agegroups.rds")
rmarkdown::paged_table(sir_data) # annotate so that it is suitable for the particle filter to use

## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I updated the code, filled the parameters with numbers;
# e.g.dt <- user(0) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_ageGroup2.R")

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

# This is part of sir odin model:
pars <- list(m = t_norm,
             N_ini = contact_2_demographic$demography$population,
             log_A_ini = 0, # c(0, 0),
             phi = 0,
             time_shift_1 = 0,
             beta_0 = 0,
             beta_1 = 0,
             # beta_diff = 0,
             vacc = 0,
             log_delta1 = 0,
             log_delta2 = 0,
             # sigma_1 = 0,
             omega = 0
)

# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model-over-time
# n_particles <- 50 # Trial n_particles = 50
# filter <- mcstate::particle_filter$new(data = sir_data,
#                                        model = gen_sir, # Use odin.dust input
#                                        n_particles = n_particles,
#                                        compare = case_compare,
#                                        seed = 1L)
# 
# filter$run(pars)

# Variance and particles estimation (as suggested by Rich)
# parallel::detectCores() # 4 cores
# x <- replicate(30, filter$run(pars))
# var(x)
# [1] 3520.937
# Trial 320000 particles to get var(x) = 1 on 4 chains/4 nodes/4 cores per-node
# (320000/4/4/4)/3520.937
# Trial 320000 particles to get var(x) = 1 on 4 chains/1 nodes/20 cores per-node
# 3520.937*(4*1*20) # ~ 281675
# (281675/4/1/20)/3520.937

# Update n_particles based on calculation in 4 cores with var(x) ~ 3520.937: 281675

priors <- prepare_priors(pars)
proposal_matrix <- diag(0.1, 9)
# diagonal ≈ (reasonable_range/k)^2
# (k = 3: extreme jump, 5 or 6 to be more conservative)
# k <- 5
# proposal_matrix <- diag(c(
#   (0.1/k)^2, # log_A_ini
#   (0.1/k)^2, # phi
#   (0.1/k)^2, # time_shift_1
#   (0.005/k)^2, # beta_0
#   (0.1/k)^2, # beta_1
#   (0.1/k)^2, # log_delta1
#   (0.1/k)^2 # log_delta2
#   # (2/5)^2 # kappa_1
# ))
rownames(proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")
colnames(proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")

mcmc_pars <- prepare_parameters(initial_pars = pars,
                                priors = priors,
                                proposal = proposal_matrix,
                                transform = transform)

# Check the transform function:
# transform(mcmc_pars$initial())

# n_steps <- 100 #1e6

# I change pmcmc_run into a function that involve control inside:
# pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)

# Directory for saving the outputs
# dir.create("outputs/genomics/trial_deterministic_1000", FALSE, TRUE)

# Trial combine pMCMC + tuning #################################################
pmcmc_run_plus_tuning <- function(n_pars, n_sts,
                                  run1_stochastic = TRUE, run2_stochastic = TRUE,
                                  ncpus){
  # dir_name <- paste0("outputs/genomics/trial_", ifelse(run_stochastic, "stochastic", "deterministic"), "_", n_sts, "/")
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(dir_name, FALSE, TRUE)
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  
  if(run1_stochastic){
    filter <- mcstate::particle_filter$new(data = sir_data,
                                           model = gen_sir, # Use odin.dust input
                                           n_particles = n_pars,
                                           compare = case_compare,
                                           seed = 1L)
  } else {
    filter <- mcstate::particle_deterministic$new(data = sir_data,
                                                  model = gen_sir,
                                                  compare = case_compare,
                                                  index = index_fun)
  }
  
  adaptive_proposal_run1 <- mcstate::adaptive_proposal_control(initial_vcv_weight = 100,
                                                               initial_scaling = (2.38^2/8), #/1e4,
                                                               # scaling_increment = NULL,
                                                               acceptance_target = 0.3,
                                                               forget_rate = 0.2,
                                                               forget_end = Inf,
                                                               adapt_end = Inf,
                                                               pre_diminish = Inf
  )
  
  # adjust short run 20% of MCMC2
  control <- mcstate::pmcmc_control(n_steps = n_sts/5,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE,
                                    
                                    n_chains = 1,
                                    # n_workers = 4,
                                    n_threads_total = ncpus,
                                    save_state = TRUE,
                                    save_trajectories = TRUE
                                    # adaptive_proposal = adaptive_proposal_run1
                                    )
  
  # The pmcmc
  pmcmc_result <- mcstate::pmcmc(mcmc_pars, filter, control = control)
  # pmcmc_result
  # saveRDS(pmcmc_result, paste0(dir_name, "pmcmc_result.rds"))
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, paste0(dir_name, "new_proposal_mtx.csv"), row.names = FALSE)
  
  lpost_max <- which.max(pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(pmcmc_result$pars[lpost_max, ]),
            paste0(dir_name, "initial.csv"), row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc1 <- pmcmc_further_process(n_sts/10, pmcmc_result)
  write.csv(mcmc1, paste0(dir_name, "mcmc1.csv"), row.names = FALSE)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc1)
  write.csv(calc_ess, paste0(dir_name, "calc_ess.csv")) #, row.names = FALSE)
  
  # Figures! (still failed, margin error)
  # # fig <- pmcmc_trace(mcmc1)
  
  Sys.sleep(1)
  
  # New proposal matrix
  new_proposal_matrix <- as.matrix(read.csv(paste0(dir_name, "new_proposal_mtx.csv")))
  new_proposal_matrix <- apply(new_proposal_matrix, 2, as.numeric)
  # vcv positive definite error if matrix/1000
  # new_proposal_matrix <- new_proposal_matrix*1.2
  new_proposal_matrix <- (new_proposal_matrix + t(new_proposal_matrix))/2
  rownames(new_proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")
  colnames(new_proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")
  # isSymmetric(new_proposal_matrix)
  
  tune_mcmc_pars <- prepare_parameters(initial_pars = pars,
                                       priors = priors,
                                       proposal = new_proposal_matrix,
                                       transform = transform)
  
  # Including adaptive proposal control
  # https://mrc-ide.github.io/mcstate/reference/adaptive_proposal_control.html
  # note:
  # use (MCMC1 & turn off adaptive_proposal) OR (just MCMC2 with adaptive_proposal)
  adaptive_proposal_run2 <- mcstate::adaptive_proposal_control(initial_vcv_weight = 1,
                                                               initial_scaling = (2.38^2/nrow(new_proposal_matrix)),
                                                               # scaling_increment = NULL,
                                                               acceptance_target = 0.234,
                                                               forget_rate = 0.2,
                                                               forget_end = Inf,
                                                               adapt_end = Inf,
                                                               pre_diminish = 0.5
  )
  
  
  if(run2_stochastic){
    tune_control <- mcstate::pmcmc_control(n_steps = n_sts,
                                           rerun_every = 50,
                                           rerun_random = TRUE,
                                           progress = TRUE,
                                           
                                           n_chains = 4,
                                           # n_workers = 4,
                                           n_threads_total = ncpus,
                                           save_state = TRUE,
                                           save_trajectories = TRUE)
    
    
    filter <- mcstate::particle_filter$new(data = sir_data,
                                           model = gen_sir, # Use odin.dust input
                                           n_particles = n_pars,
                                           compare = case_compare,
                                           seed = 1L
    )
  } else {
    tune_control <- mcstate::pmcmc_control(n_steps = n_sts,
                                           rerun_every = 50,
                                           rerun_random = TRUE,
                                           progress = TRUE,
                                           
                                           n_chains = 4,
                                           # n_workers = 4,
                                           n_threads_total = ncpus,
                                           save_state = TRUE,
                                           save_trajectories = TRUE
                                           # another option is to construct vcv first then ignore adaptive_proposal settings
                                           # adaptive_proposal = adaptive_proposal_run2
    )
    
    filter <- mcstate::particle_deterministic$new(data = sir_data,
                                                  model = gen_sir, # Use odin.dust input
                                                  compare = case_compare,
                                                  index = index_fun
    )
  }
  
  # The pmcmc
  tune_pmcmc_result <- mcstate::pmcmc(tune_mcmc_pars, filter, control = tune_control)
  # saveRDS(tune_pmcmc_result, paste0(dir_name, "tune_pmcmc_result.rds"))
  
  # final parameters with CI
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  mcmc_lo_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.025))
  mcmc_hi_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.975))
  
  binds_tune_initial <- rbind(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
                              mcmc_lo_CI, mcmc_hi_CI)
  binds_tune_initial2 <- cbind(binds_tune_initial,
                               log_prior = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                           "log_prior"],
                               log_likelihood = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                                "log_likelihood"],
                               log_posterior = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                               "log_posterior"])
  t_tune_initial <- t(binds_tune_initial2)
  colnames(t_tune_initial) <- c("values", "low_CI", "high_CI")
  
  write.csv(t_tune_initial,
            paste0(dir_name, "tune_initial_with_CI.csv"), row.names = T)
  
  # MCMC diagnostics
  # 1. Gelman-Rubin
  figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  
  png(paste0(dir_name, "figs/mcmc2_diag_gelmanRubin_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  diag_gelman_rubin(figs_gelman_init)
  dev.off()
  
  # 2. ggpairs
  png(paste0(dir_name, "figs/mcmc2_diag_ggPairs_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  p <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  print(p)
  dev.off()
  
  new_proposal_mtx <- cov(tune_pmcmc_result$pars)
  write.csv(new_proposal_mtx, paste0(dir_name, "new_proposal_mtx_modified.csv"), row.names = FALSE)
  
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
            paste0(dir_name, "tune_initial.csv"), row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc2 <- coda::as.mcmc(cbind(
    tune_pmcmc_result$probabilities, tune_pmcmc_result$pars))
  write.csv(mcmc2, paste0(dir_name, "mcmc2.csv"), row.names = FALSE)
  
  mcmc2_burnedin <- tuning_pmcmc_further_process(n_sts, tune_pmcmc_result)
  write.csv(mcmc2_burnedin, paste0(dir_name, "mcmc2_burnedin.csv"), row.names = FALSE)
  
  # Calculating ESS & Acceptance Rate
  tune_calc_ess <- ess_calculation(mcmc2)
  write.csv(tune_calc_ess, paste0(dir_name, "tune_calc_ess.csv")) #, row.names = FALSE)
  
  tune_calc_ess_burnedin <- ess_calculation(mcmc2_burnedin)
  write.csv(tune_calc_ess_burnedin, paste0(dir_name, "tune_calc_ess_burnedin.csv")) #, row.names = FALSE)
  
  # Figures! (still failed, margin error)
  # fig <- pmcmc_trace(mcmc2)
  # fig <- pmcmc_trace(mcmc2_burnedin)
  
  # save pmcmc samples
  if(run2_stochastic){
    pmcmc_samples <- mcstate::pmcmc_sample(tune_pmcmc_result,
                                           n_sample = 1000,
                                           burnin = 500)
    
    saveRDS(pmcmc_samples, paste0(dir_name, "pmcmc_samples.rds"))
  } else {
    # message("Not required for a deterministic model")
    pmcmc_samples <- mcstate::pmcmc_sample(tune_pmcmc_result,
                                           n_sample = 1000,
                                           burnin = 500)
    
    saveRDS(pmcmc_samples, paste0(dir_name, "pmcmc_samples.rds"))
  }
  
  
  ##############################################################################
  # MCMC Diagnostics
  
  # 1. Gelman-Rubin Diagnostic
  # https://cran.r-project.org/web/packages/coda/coda.pdf
  figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  # fig <- diag_cov_mtx(figs_gelman_init)
  # fig <- diag_gelman_rubin(figs_gelman_init)
  
  # 2. Autocorrelation
  # fig <- diag_aucorr(mcmc2)
  
  # 3. ggpairs
  # fig <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  
}


# run pmcmc1 only
pmcmc_run1_only <- function(n_pars, n_sts,
                            run1_stochastic = TRUE,
                            ncpus){
  # dir_name <- paste0("outputs/genomics/trial_", ifelse(run_stochastic, "stochastic", "deterministic"), "_", n_sts, "/")
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(dir_name, FALSE, TRUE)
  
  if(run1_stochastic){
    filter <- mcstate::particle_filter$new(data = sir_data,
                                           model = gen_sir, # Use odin.dust input
                                           n_particles = n_pars,
                                           compare = case_compare,
                                           seed = 1L)
  } else {
    filter <- mcstate::particle_deterministic$new(data = sir_data,
                                                  model = gen_sir,
                                                  compare = case_compare,
                                                  index = index_fun)
  }
  
  
  
  control <- mcstate::pmcmc_control(n_steps = n_sts,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE,
                                    
                                    n_chains = 1,
                                    # n_workers = 4,
                                    n_threads_total = ncpus,
                                    save_state = TRUE,
                                    save_trajectories = TRUE,
                                    adaptive_proposal = adaptive_proposal_run1
                                    )
  
  # The pmcmc
  pmcmc_result <- mcstate::pmcmc(mcmc_pars, filter, control = control)
  pmcmc_result
  saveRDS(pmcmc_result, paste0(dir_name, "pmcmc_result.rds"))
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, paste0(dir_name, "new_proposal_mtx.csv"), row.names = FALSE)
  
  lpost_max <- which.max(pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(pmcmc_result$pars[lpost_max, ]),
            paste0(dir_name, "initial.csv"), row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc1 <- pmcmc_further_process(n_sts, pmcmc_result)
  write.csv(mcmc1, paste0(dir_name, "mcmc1.csv"), row.names = FALSE)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc1)
  write.csv(calc_ess, paste0(dir_name, "calc_ess.csv")) #, row.names = FALSE)
  
  # Figures! (still failed, margin error)
  # # fig <- pmcmc_trace(mcmc1)
  
}


# run pmcmc2 only
pmcmc_run2_only <- function(n_pars, n_sts,
                             run2_stochastic = TRUE,
                             ncpus){
  # dir_name <- paste0("outputs/genomics/trial_", ifelse(run_stochastic, "stochastic", "deterministic"), "_", n_sts, "/")
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(dir_name, FALSE, TRUE)
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  
  # New proposal matrix is self-designed diagonal matrix
  # proposal_matrix <- diag(0.1, 7)
  # diagonal ≈ (reasonable_range/k)^2
  # (k = 3: extreme jump, 5 or 6 to be more conservative)
  k <- 5
  proposal_matrix <- diag(c(
    (0.08/k)^2, # log_A_ini
    (0.25/k)^2, # phi
    (0.08/k)^2, # time_shift_1
    (0.002/k)^2, # beta_0; quite sensitive must be < 0.005
    (0.12/k)^2, # beta_1
    # (0.1/k)^2, # beta_diff
    (3e-5/k)^2, # vacc
    (0.08/k)^2, # log_delta1
    (0.08/k)^2, # log_delta2
    # (0.01/k)^2, # sigma_1
    (2e-4/k)^2 # omega
    # (2/k)^2 # kappa_1
  ))
  
  new_proposal_matrix <- as.matrix(proposal_matrix)
  new_proposal_matrix <- (new_proposal_matrix + t(new_proposal_matrix))/2
  rownames(new_proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")
  colnames(new_proposal_matrix) <- c("log_A_ini", "phi", "time_shift_1", "beta_0", "beta_1", "vacc", "log_delta1", "log_delta2", "omega")
  # isSymmetric(new_proposal_matrix)
  
  tune_mcmc_pars <- prepare_parameters(initial_pars = pars,
                                       priors = priors,
                                       proposal = new_proposal_matrix,
                                       transform = transform)
  
  # Including adaptive proposal control
  # https://mrc-ide.github.io/mcstate/reference/adaptive_proposal_control.html
  # note:
  # use (MCMC1 & turn off adaptive_proposal) OR (just MCMC2 with adaptive_proposal)
  adaptive_proposal_run2 <- mcstate::adaptive_proposal_control(initial_vcv_weight = 5, # 1 for non-vaccine model
                                                               initial_scaling = (2.38^2/nrow(new_proposal_matrix))*0.5,
                                                               scaling_increment = 0.05,
                                                               acceptance_target = 0.234,
                                                               forget_rate = 0.2,
                                                               forget_end = Inf,
                                                               adapt_end = Inf,
                                                               pre_diminish = 0.5
  )
  
  
  if(run2_stochastic){
    tune_control <- mcstate::pmcmc_control(n_steps = n_sts,
                                           rerun_every = 50,
                                           rerun_random = TRUE,
                                           progress = TRUE,
                                           
                                           n_chains = 4,
                                           # n_workers = 4,
                                           n_threads_total = ncpus,
                                           save_state = TRUE,
                                           save_trajectories = TRUE)
    
    
    filter <- mcstate::particle_filter$new(data = sir_data,
                                           model = gen_sir, # Use odin.dust input
                                           n_particles = n_pars,
                                           compare = case_compare,
                                           seed = 1L
    )
  } else {
    tune_control <- mcstate::pmcmc_control(n_steps = n_sts,
                                           rerun_every = 50,
                                           rerun_random = TRUE,
                                           progress = TRUE,
                                           
                                           n_chains = 4,
                                           # n_workers = 4,
                                           n_threads_total = ncpus,
                                           save_state = TRUE,
                                           save_trajectories = TRUE,
                                           # another option is to construct vcv first then ignore adaptive_proposal settings
                                           adaptive_proposal = adaptive_proposal_run2
    )
    
    filter <- mcstate::particle_deterministic$new(data = sir_data,
                                                  model = gen_sir, # Use odin.dust input
                                                  compare = case_compare,
                                                  index = index_fun
    )
  }
  
  # The pmcmc
  tune_pmcmc_result <- mcstate::pmcmc(tune_mcmc_pars, filter, control = tune_control)
  # saveRDS(tune_pmcmc_result, paste0(dir_name, "tune_pmcmc_result.rds"))
  
  # final parameters with CI
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  mcmc_lo_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.025))
  mcmc_hi_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.975))
  
  binds_tune_initial <- rbind(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
                              mcmc_lo_CI, mcmc_hi_CI)
  binds_tune_initial2 <- cbind(binds_tune_initial,
                               log_prior = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                           "log_prior"],
                               log_likelihood = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                                "log_likelihood"],
                               log_posterior = tune_pmcmc_result$probabilities[tune_lpost_max,
                                                                               "log_posterior"])
  t_tune_initial <- t(binds_tune_initial2)
  colnames(t_tune_initial) <- c("values", "low_CI", "high_CI")
  
  write.csv(t_tune_initial,
            paste0(dir_name, "tune_initial_with_CI.csv"), row.names = T)
  
  # MCMC diagnostics
  # 1. Gelman-Rubin
  figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  
  png(paste0(dir_name, "figs/mcmc2_diag_gelmanRubin_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  diag_gelman_rubin(figs_gelman_init)
  dev.off()
  
  # 2. ggpairs
  png(paste0(dir_name, "figs/mcmc2_diag_ggPairs_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  p <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  print(p)
  dev.off()
  
  new_proposal_mtx <- cov(tune_pmcmc_result$pars)
  write.csv(new_proposal_mtx, paste0(dir_name, "new_proposal_mtx_modified.csv"), row.names = FALSE)
  
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]),
            paste0(dir_name, "tune_initial.csv"), row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc2 <- coda::as.mcmc(cbind(
    tune_pmcmc_result$probabilities, tune_pmcmc_result$pars))
  write.csv(mcmc2, paste0(dir_name, "mcmc2.csv"), row.names = FALSE)
  
  mcmc2_burnedin <- tuning_pmcmc_further_process(n_sts, tune_pmcmc_result)
  write.csv(mcmc2_burnedin, paste0(dir_name, "mcmc2_burnedin.csv"), row.names = FALSE)
  
  # Calculating ESS & Acceptance Rate
  tune_calc_ess <- ess_calculation(mcmc2)
  write.csv(tune_calc_ess, paste0(dir_name, "tune_calc_ess.csv")) #, row.names = FALSE)
  
  tune_calc_ess_burnedin <- ess_calculation(mcmc2_burnedin)
  write.csv(tune_calc_ess_burnedin, paste0(dir_name, "tune_calc_ess_burnedin.csv")) #, row.names = FALSE)
  
  # Figures! (still failed, margin error)
  # fig <- pmcmc_trace(mcmc2)
  # fig <- pmcmc_trace(mcmc2_burnedin)
  
  # save pmcmc samples
  if(run2_stochastic){
    pmcmc_samples <- mcstate::pmcmc_sample(tune_pmcmc_result,
                                           n_sample = 1000,
                                           burnin = 500)
    
    saveRDS(pmcmc_samples, paste0(dir_name, "pmcmc_samples.rds"))
  } else {
    # message("Not required for a deterministic model")
    pmcmc_samples <- mcstate::pmcmc_sample(tune_pmcmc_result,
                                           n_sample = 1000,
                                           burnin = 500)
    
    saveRDS(pmcmc_samples, paste0(dir_name, "pmcmc_samples.rds"))
  }
  
}

# a slight modification for college's HPC
# PRECAUTION: will be resulted in errors if run locally
args <- commandArgs(trailingOnly = T)

# get_arg as a function
get_arg <- function(flag){
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    stop(paste("Missing/malformed argument for", flag))
  }
  args[idx + 1]
}

mode <- get_arg("--mode")
n_pars <- as.numeric(get_arg("--n_particles"))
n_sts <- as.numeric(get_arg("--n_steps"))
ncpus <- as.numeric(get_arg("--ncpus"))

run1_stochastic_flag <- tolower(get_arg("--run1_stochastic"))
run1_stochastic <- run1_stochastic_flag %in% c("true", "t", "1")
if (mode %in% c("pmcmc2", "run_all")) {
  run2_stochastic_flag <- tolower(get_arg("--run2_stochastic"))
  run2_stochastic <- run2_stochastic_flag %in% c("true", "t", "1")
} else {
  run2_stochastic <- NA
}


if (mode == "pmcmc1") {
  pmcmc_run1_only(n_pars, n_sts, run1_stochastic, ncpus)
} else if (mode == "pmcmc2") {
  pmcmc_run2_only(n_pars, n_sts, run2_stochastic, ncpus)
} else if (mode == "run_all") {
  pmcmc_run_plus_tuning(n_pars, n_sts, run1_stochastic, run2_stochastic, ncpus)
} else {
  stop("Invalid mode selected. Choose: pmcmc1, pmcmc2, or run_all.")
}


# pmcmc_run1_only(10, 600, run1_stochastic = T, ncpus = 4)
# pmcmc_run2_only(10, 600, run2_stochastic = F, ncpus = 4)
# pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
#                       run1_stochastic = F, run2_stochastic = F, ncpus = 4)

# best results: turn on rerun_every & rerun_random during pmcmc2 but do not use adaptive proposal control
