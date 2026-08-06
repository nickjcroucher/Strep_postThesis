# analyse pmcmc samples
# related fun

source("global/all_function_allAge.R")

post_particle_pics <- function(n_sts){
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  
  pmcmc_samples <- readRDS(paste0(dir_name, "pmcmc_samples.rds"))
  pmcmc_samples$trajectories$state <- observe(pmcmc_samples)
  
  data <- readRDS("raw_data/pmcmc_data_week_allAge_ser1_test_2agegroups.rds")
  
  if(file.exists(file.path(paste0(dir_name, "tune_initial.csv")))){
    initial_pars <- read.csv(paste0(dir_name, "tune_initial.csv"))
    priors <- prepare_priors(initial_pars)
    
    png(paste0(dir_name, "figs/particles_posteriors_%02d.png"),
        width = 24, height = 17, unit = "cm", res = 600)
    par(mfrow = c(3, 3), mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0), bty = "n")
    for (nm in names(initial_pars)) {
      hist(pmcmc_samples$pars[, nm], xlab = nm, main = "", freq = FALSE)
      if (nm %in% names(priors)) {
        curve(exp(priors[[nm]](x)), add = TRUE, col = 2)
      }
    }
    dev.off()
  }
  
  png(paste0(dir_name, "figs/particles_fits_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  par(mfrow = c(2, 1), bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0))
  plot_states(pmcmc_samples$trajectories$state, data)
  dev.off()
  
  png(paste0(dir_name, "figs/mcmc2_diag_ggPairs_particles_%02d.png"),
      width = 17, height = 17, unit = "cm", res = 600)
  par(mar = c(3, 3, 1, 1), mgp = c(1.7, 0.7, 0))
  p <- GGally::ggpairs(as.data.frame(pmcmc_samples$pars))
  print(p)
  dev.off()
  
}


# a slight modification for college's HPC
args <- commandArgs(trailingOnly = T)
n_sts <- as.numeric(args[which(args == "--n_steps") + 1])

post_particle_pics(n_sts)
