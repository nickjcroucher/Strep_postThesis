rm(list = ls())
library(tidyverse)

if (!dir.exists("inputs")) {
  dir.create("inputs")
}

source("global/all_function_allAge.R")
# global/all_function_allAge.R also incorporated:
# burnin_days

# Data preparation for serotype 1
# mcstate data preparation #####################################################
# non-heterogeneity (allAges), weekly
serotype1_data <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv") %>% 
  glimpse()

all_week <- data.frame(week_date = seq.Date(from = min(as.Date(serotype1_data$Earliest.specimen.date)),
                                            to = max(as.Date(serotype1_data$Earliest.specimen.date)), 
                                            by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
                ) %>% 
  distinct(yearWeek) %>% 
  glimpse()


allAges_weekly_ser1 <- dplyr::left_join(
  all_week
  ,
  serotype1_data %>% 
    dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                  iso_week = paste0(year(Earliest.specimen.date), "-W", sprintf("%02d", week(Earliest.specimen.date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week),
                  
                  ageGroup_s1 = case_when(
                    AGEYR < 15 ~ "0-14",
                    AGEYR >= 15 ~ "15+"
                  )
    ) %>% 
    dplyr::filter(!is.na(ageGroup_s1)) %>% 
    dplyr::group_by(yearWeek, ageGroup_s1) %>% 
    dplyr::summarise(count_serotype = n()) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(
      .,
      names_from = contains("ageGroup"),
      names_prefix = "count_",
      values_from = "count_serotype"
    ) %>% 
    dplyr::rename(
      count_s1_1 = "count_0-14",
      count_s1_2 = "count_15+"
    ) %>% 
    dplyr::arrange(yearWeek)
  ,
  by = "yearWeek"
  ) %>% 
  dplyr::mutate(
    yearWeek = as.Date(yearWeek),
    day = as.numeric(round((yearWeek - (as.Date("2003-01-01")-2)))), # min(dat_g$Earliest.specimen.date)-2 to make it 7
    
    # fill in NAs
    # count_s1_1 = ifelse(is.na(count_s1_1), 0, count_s1_1),
    # count_s1_2 = ifelse(is.na(count_s1_2), 0, count_s1_2),
    
    # adjust burn in
    day = burnin_days+day,
    count_s1_1 = ifelse(dplyr::row_number() == 1, NA, count_s1_1),
    count_s1_2 = ifelse(dplyr::row_number() == 1, NA, count_s1_2)
  ) %>% 
  dplyr::filter(day > 0,
                # adjust data fitting from the lowest case peak post-winter
                # yearWeek >= as.Date("2003-03-01")
                ) %>%
  mcstate::particle_filter_data(.,
                                time = "day", # I use steps instead of day
                                rate = 1, # I change the model to weekly, therefore weekly rate is required
                                initial_time = 0
  ) %>%
  glimpse()

saveRDS(allAges_weekly_ser1, "raw_data/pmcmc_data_week_allAge_ser1_test_2agegroups.rds")

# test plot
ggplot(allAges_weekly_ser1
       , aes(x = yearWeek)) +
  geom_line(size = 0.5, aes(y = count_s1_1), colour = "darkgreen") +
  geom_line(size = 0.5, aes(y = count_s1_2), colour = "maroon") +
  geom_vline(xintercept = as.Date("2010-04-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2002-12-31"), as.Date("2020-12-31")),
               date_breaks = "1 year",
               date_labels = "%Y") +
  # scale_y_log10() +
  theme_bw() +
  labs(
    title = "Serotype 1 Counts",
    y = "Serotype 1 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))










