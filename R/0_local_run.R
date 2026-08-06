source("global/all_function_allAge.R")
source("model/sir_stochastic_ageGroup2_2post_pmcmc_picts.R")

source("R/3_pmcmc.R")
source("R/4_post_pmcmc_pics.R")
source("R/5_post_pmcmc_samples_pics.R")
source("R/6_post_pmcmc_age_validation.R")

# test vcv error
count <- 4
repeat { 
  pmcmc_run_plus_tuning(n_pars = 10, n_sts = 100,
                        run1_stochastic = F, run2_stochastic = F, ncpus = 60)
  if (count < 5)
    break
}

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)
post_pmcmc_pics(600)
model_vs_data(600)
post_particle_pics(600)
age_validation(600)


# vcv chunk error retry
max_retries <- 5

run_chunk <- function(n_sts, attempt = 1) {
  message(sprintf("Running %d (attempt %d)", n_sts, attempt))
  
  tryCatch({
    pmcmc_run_plus_tuning(n_pars = 10, n_sts = n_sts,
                          run1_stochastic = F, run2_stochastic = F, ncpus = 60)
    
    post_pmcmc_pics(n_sts)
    model_vs_data(n_sts)
    post_particle_pics(n_sts)
    age_validation(n_sts)
    
  }, error = function(e) {
    message(sprintf("Error in %d: %s", n_sts, e$message))
    
    if (attempt < max_retries) {
      message(sprintf("Retrying... (%d/%d)", attempt + 1, max_retries))
      Sys.sleep(1)
      run_chunk(n_sts, attempt + 1)
    } else {
      message(sprintf("n_sts = %d failed after %d attempts — skipping", 
                      n_sts, max_retries))
    }
  })
}


for (n in c(5000, 30000)){ #, 100000)){ 5000, 10000, 20000, 
  run_chunk(n)
}
