library(tidyverse)
library(odin.dust)
library(socialmixr)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")
source("global/all_function_allAge.R")
# global/all_function_allAge.R also incorporated:
# burnin_days

vaccine_simulation <- function(vaccYear){
  dir_name <- paste0("outputs/genomics/trial_", 500000, "/")
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  # run 4_post_pmcmc_pics.R first
  results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv"),
                      row.names = 1) %>% 
    glimpse()
  
  # adjust year <-> days based on vaccYear
  vdays <- switch(
    as.character(vaccYear),
    "2003" = 0,
    "2006" = 1339,
    "2010" = 2648,
    "2023" = 7364, # basically no vaccine being introduced
    stop("Invalid year")
  )
  
  # gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
  gen_sir <- odin.dust::odin_dust("model/sir_stochastic_ageGroup2_adjusted_vacc.R")
  
  # Create contact_matrix 5 demographic groups:
  # > 5
  # 5-18
  # 19-30
  # 31-64
  # 65+
  # age.limits = c(0, 5, 19, 31, 65)
  
  # Create contact_matrix 2 demographic groups:
  # < 15
  # 15+
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
  
  pars <- list(m = t_norm,
               N_ini = contact_2_demographic$demography$population,
               vaccIntro = as.numeric(vdays),
               log_A_ini = results["log_A_ini",2], # c(results[1,2], results[2,2]),
               phi = results["phi",2],
               time_shift_1 = results["time_shift_1",2],
               beta_0 = results["beta_0",2],
               beta_1 = results["beta_1",2],
               vacc = results["vacc",2],
               log_delta1 = results["log_delta1",2],
               # rho = results["rho",2],
               log_delta2 = results["log_delta2",2],
               # sigma_1 = results["sigma_1",2],
               omega = results["omega",2]
  )
  
  n_times <- burnin_days+7500
  n_pars <- 10000L
  sir_model <- gen_sir$new(pars = pars,
                           time = 1,
                           n_particles = n_pars,
                           n_threads = 4L,
                           seed = 1L)
  
  model <- array(NA, dim = c(sir_model$info()$len, n_pars, n_times))
  
  for (t in seq_len(n_times)) {
    model[ , , t] <- sir_model$run(t)
  }
  
  data <- readRDS("raw_data/pmcmc_data_week_allAge_ser1_test_2agegroups.rds") %>% 
    glimpse()
  
  sir_data <- dplyr::bind_rows(
    data %>% 
      dplyr::transmute(
        replicate = 1,
        # steps = time_start+1,
        weekly = seq_along(replicate),
        value = count_s1_1,
        compartment = "data_count_s1_1"
      )
    ,
    data %>% 
      dplyr::transmute(
        replicate = 1,
        # steps = time_start+1,
        weekly = seq_along(replicate),
        value = count_s1_2,
        compartment = "data_count_s1_2"
      )
  ) %>%
    tidyr::complete(weekly, compartment,
                    fill = list(value = 0)) %>% 
    glimpse()
  
  all_dates <- data %>%
    dplyr::select(yearWeek) %>% 
    dplyr::mutate(
      weekly = seq_along(yearWeek)
    ) %>%
    glimpse()
  
  # focused on n_AD_weekly (already in weeks)
  incidence_modelled <- 
    reshape2::melt(model) %>% 
    dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                  replicate = Var2, # Var2 = particles
                  steps = Var3       # Var3 = steps are in days, but n_AD_weekly is aggregated in weeks
    ) %>% 
    # adjust burn in
    dplyr::filter(steps > burnin_days) %>% 
    dplyr::mutate(steps = steps-burnin_days) %>% 
    # dplyr::filter(index > 8) %>%
    dplyr::mutate(compartment = 
                    dplyr::case_when(index == 1 ~ "Time",
                                     index == 2 ~ "total N",
                                     index == 3 ~ "total S",
                                     index == 4 ~ "total A",
                                     index == 5 ~ "total D",
                                     index == 6 ~ "total R",
                                     index == 7 ~ "n_AD1_weekly",
                                     index == 8 ~ "n_AD2_weekly",
                                     
                                     index == 9 ~ "model_S1",
                                     index == 10 ~ "model_S2",
                                     index == 11 ~ "A <15",
                                     index == 12 ~ "A 15+",
                                     index == 13 ~ "model_D1",
                                     index == 14 ~ "model_D2",
                                     index == 15 ~ "R <15",
                                     index == 16 ~ "R 15+"
                    )) %>% 
    dplyr::filter(index > 8) %>%
    dplyr::select(-index) %>%
    dplyr::mutate(weekly = ceiling((steps-1)/7)) %>% 
    dplyr::group_by(replicate, weekly, compartment) %>% 
    dplyr::summarise(value = sum(value, na.rm = T),
                     # date = max(date),
                     .groups = "drop") %>% 
    dplyr::ungroup() %>% 
    dplyr::bind_rows(sir_data) %>%
    tidyr::complete(weekly,
                    fill = list(value = 0)
    ) %>%
    dplyr::full_join(
      all_dates
      ,
      by = "weekly"
    ) %>%
    dplyr::filter(!is.na(yearWeek)) # %>%
  # glimpse()
  
  # group_by weekly, compartment and calculate median for ALL simulated particles
  incidence_modelled_med <- incidence_modelled %>% 
    dplyr::filter(
      compartment %in% c("model_D1")
    ) %>% 
    dplyr::group_by(yearWeek) %>% 
    dplyr::summarise(
      model_D1_med = median(value),
      model_D1_lo  = quantile(value, 0.025),
      model_D1_up  = quantile(value, 0.975),
      
      .groups = "drop"
    ) %>% 
    dplyr::full_join(
      incidence_modelled %>% 
        dplyr::filter(
          compartment %in% c("model_D2")
        ) %>% 
        dplyr::group_by(yearWeek) %>% 
        dplyr::summarise(
          model_D2_med = median(value),
          model_D2_lo  = quantile(value, 0.025),
          model_D2_up  = quantile(value, 0.975),
          
          .groups = "drop"
        )
      ,
      by = "yearWeek"
    ) %>% 
    # combine with S
    dplyr::full_join(
      incidence_modelled %>% 
        dplyr::filter(
          compartment %in% c("model_S1")
        ) %>% 
        dplyr::group_by(yearWeek) %>% 
        dplyr::summarise(
          model_S1_med = median(value),
          model_S1_lo  = quantile(value, 0.025),
          model_S1_up  = quantile(value, 0.975),
          
          .groups = "drop"
        )
      ,
      by = "yearWeek"
    ) %>% 
    dplyr::full_join(
      incidence_modelled %>% 
        dplyr::filter(
          compartment %in% c("model_S2")
        ) %>% 
        dplyr::group_by(yearWeek) %>% 
        dplyr::summarise(
          model_S2_med = median(value),
          model_S2_lo  = quantile(value, 0.025),
          model_S2_up  = quantile(value, 0.975),
          
          .groups = "drop"
        )
      ,
      by = "yearWeek"
    ) %>% 
    tidyr::pivot_longer(
      cols = contains("model_"),
      names_to = "compartment",
      values_to = "value"
    ) %>% 
    glimpse()
  
  # combine median back to incidence_modelled
  incidence_modelled2 <- dplyr::bind_rows(
    incidence_modelled,
    incidence_modelled_med
  ) %>% 
    tidyr::pivot_wider(
      id_cols = c("yearWeek", "replicate"),
      names_from = compartment,
      values_from = value,
      # values_fill = 0
    ) %>% 
    glimpse()
  
  # save ONLY the med df
  write.csv(incidence_modelled_med,
            paste0(dir_name, "incidence_modelled_serotype1_median_simulated_", vaccYear, ".csv"),
            row.names = FALSE)
  
  png(paste0(dir_name, "figs/model_vs_data_ageGroups1_median_simulated_", vaccYear, ".png"),
      width = 24, height = 14, unit = "cm", res = 600)
  p1 <- ggplot(incidence_modelled2
               ,
               aes(x = yearWeek,
                   group = interaction(replicate))) +
    geom_line(aes(y = model_D1,
                  colour = "Simulated results"),
              linewidth = 0.01,
              alpha = 0.1
    ) +
    # geom_ribbon(
    #   aes(ymin = model_D1_lo, ymax = model_D1_up),
    #   # fill = "steelblue",
    #   alpha = 0.01
    # ) +
    geom_line(aes(y = data_count_s1_1,
                  colour = "Cases")
    ) +
    geom_line(aes(y = model_D1_med,
                  colour = "Median"),
              linewidth = 0.75
    ) +
    scale_y_continuous(
      limits = c(0, 20)
    ) +
    scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
                 date_breaks = "year",
                 date_labels = "%Y") +
    scale_colour_manual(
      values = c("Simulated results" = "grey40",
                 "Median" = "black",
                 "Cases" = "#FF7F00"
      )
    ) +
    ggtitle(paste0("Simulated cases for vaccine introduction in ", vaccYear, " in age group 0-14")) +
    xlab("Time") +
    ylab("Number of People") +
    theme_bw() +
    theme(legend.position = c(0.15, 0.85),
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "lines"),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "transparent", color = "transparent"))
  
  p2 <- ggplot(incidence_modelled2
               ,
               aes(x = yearWeek,
                   group = interaction(replicate))) +
    geom_line(aes(y = model_D2,
                  colour = "Simulated results"),
              linewidth = 0.01,
              alpha = 0.1
    ) +
    # geom_ribbon(
    #   aes(ymin = model_D2_lo, ymax = model_D2_up),
    #   # fill = "steelblue",
    #   alpha = 0.01
    # ) +
    geom_line(aes(y = data_count_s1_2,
                  colour = "Cases")
    ) +
    geom_line(aes(y = model_D2_med,
                  colour = "Median"),
              linewidth = 0.75
    ) +
    scale_y_continuous(
      limits = c(0, 50)
    ) +
    scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
                 date_breaks = "year",
                 date_labels = "%Y") +
    scale_colour_manual(
      values = c("Simulated results" = "grey40",
                 "Median" = "black",
                 "Cases" = "#FF7F00"
      )
    ) +
    ggtitle(paste0("Simulated cases for vaccine introduction in ", vaccYear, " in age group 15+")) +
    xlab("Time") +
    ylab("Number of People") +
    theme_bw() +
    theme(legend.position = c(0.15, 0.85),
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "lines"),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "transparent", color = "transparent"))
  
  p_combined <- cowplot::plot_grid(p1, p2,
                                   nrow =2,
                                   labels = c("A", "B"))
  
  
  print(p_combined)
  dev.off()
  
}

# a slight modification for college's HPC
args <- commandArgs(trailingOnly = T)
vaccYear <- as.numeric(args[which(args == "--year") + 1])

vaccine_simulation(vaccYear)

