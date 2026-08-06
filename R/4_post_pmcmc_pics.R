source("global/all_function_allAge.R")

# generate post-mcmc picts
post_pmcmc_pics <- function(n_sts){
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  
  if(file.exists(file.path(paste0(dir_name, "mcmc1.csv")))){
    mcmc1 <- read.csv(paste0(dir_name, "mcmc1.csv"))
    
    # mcmc1
    # fig <- pmcmc_trace(coda::as.mcmc(mcmc1))
    png(paste0(dir_name, "figs/mcmc1_%02d.png"),
        width = 17, height = 17, unit = "cm", res = 600)
    pmcmc_trace(coda::as.mcmc(mcmc1))
    dev.off()
  }
  
  if(file.exists(file.path(paste0(dir_name, "mcmc2.csv")))){
    mcmc2 <- read.csv(paste0(dir_name, "mcmc2.csv"))
    mcmc2_burnedin <- read.csv(paste0(dir_name, "mcmc2_burnedin.csv"))
    # tune_pmcmc_result <- readRDS(paste0(dir_name, "tune_pmcmc_result.rds"))
    
    # mcmc2
    # fig <- pmcmc_trace(coda::as.mcmc(mcmc2))
    png(paste0(dir_name, "figs/mcmc2_%02d.png"),
        width = 17, height = 17, unit = "cm", res = 600)
    pmcmc_trace(coda::as.mcmc(mcmc2))
    dev.off()
    
    # mcmc2 burned in
    # fig <- pmcmc_trace(coda::as.mcmc(mcmc2))
    png(paste0(dir_name, "figs/mcmc2_burnedin_%02d.png"),
        width = 17, height = 17, unit = "cm", res = 600)
    pmcmc_trace(coda::as.mcmc(mcmc2_burnedin))
    dev.off()
    
    # 2. Autocorrelation
    png(paste0(dir_name, "figs/mcmc2_diag_auCorr_%02d.png"),
        width = 17, height = 17, unit = "cm", res = 600)
    diag_aucorr(coda::as.mcmc(mcmc2))
    dev.off()
  }
}

# a slight modification for college's HPC
args <- commandArgs(trailingOnly = T)
n_sts <- as.numeric(args[which(args == "--n_steps") + 1])

post_pmcmc_pics(n_sts)

