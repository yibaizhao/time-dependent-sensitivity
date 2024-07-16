library(here)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)

# datastamp <- "2024-04-15"
# datastamp <- "2024-04-17"
# datastamp <- "2024-04-24"
# datastamp <- "2024-04-25"
# datastamp <- "2024-05-01"
# datastamp <- "2024-05-07"
# datastamp <- "2024-05-01"
# datastamp <- "2024-07-15"
datastamp <- "2024-07-16"

##################################################
# True sensitivity evaluated at specified times
# sojourn_time: sojourn time in years
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
# TODO: Determine annotation positions programmatically
##################################################
test_sensitivity <- function(time,
                             sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             is_indolent=FALSE,
                             method='linear'){
  if(!is_indolent){
    if(method=='linear'){
      alpha <- onset_sensitivity
      beta <- (clinical_sensitivity-onset_sensitivity)/sojourn_time
      tset <- tibble(time=time, sensitivity=alpha+beta*time)
    }else{
      alpha <- -log(clinical_sensitivity/onset_sensitivity)
      beta <- -log(onset_sensitivity/clinical_sensitivity)-alpha
      tset <- tibble(time=time, sensitivity=1/(1+exp(-(alpha+beta/sojourn_time*time))))
    }
    tset <- tset %>% mutate(sensitivity = ifelse(time>=sojourn_time, clinical_sensitivity, 
                                                 ifelse(time<0, onset_sensitivity, sensitivity)))
  }else{
    tset <- tibble(time=time, 
                   sensitivity=onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/sojourn_time*time)))
  }
  return(tset)
}

# calculate P(tc<=C)
CDF_tp_tc <- function(t0, C, lambda1, lambda2){
  z <- C - t0
  cdf_Z <- 1 - (lambda2/(lambda2-lambda1))*exp(-lambda1*z) + (lambda1/(lambda2-lambda1))*exp(-lambda2*z)
  
  return(cdf_Z)
}
# calculate P(tp<=C)
CDF_tp <- function(t0, C, lambda1){
  x <- C - t0
  cdf_X <- 1 - exp(-lambda1*x)
  
  return(cdf_X)
}

# calculate tp<=K<tc, K is test age
CDF_tp_K_tc <- function(t0, K, lambda1, lambda2){
  CDF_tp_k <- CDF_tp(t0, K, lambda1)
  CDF_tc_k_condition <- lambda1*exp(-lambda2*K) * (exp((lambda2-lambda1)*K) - 1)/(lambda2-lambda1)
  
  return(CDF_tp_k * CDF_tc_k_condition)
}

# calculate t_p<=C and t_c<=C
summary_table <- function(mean_onset_time, 
                          mean_sojourn_time,
                          threshold,
                          start_age,
                          test_age){
  lambda1 <- 1/mean_onset_time # parameter for onset time ~ exp(lambda1)
  lambda2 <- 1/mean_sojourn_time # parameter for sojourn time ~ exp(lambda2)

  sset <- expand.grid(start_age = start_age,
                      test_age_K = test_age,
                      mean_onset_time = mean_onset_time, 
                      mean_sojourn_time = mean_sojourn_time,
                      threshold_C = threshold)
  
  sset <- sset %>% 
    mutate(prob_tp_C = CDF_tp(t0 = start_age, C = threshold_C, 
                              lambda1 = 1/mean_onset_time),
           prob_tc_C = CDF_tp_tc(t0 = start_age, C = threshold_C, 
                                 lambda1 = 1/mean_onset_time, lambda2 = 1/mean_sojourn_time),
           prob_tp_K_tc = CDF_tp_K_tc(t0 = start_age, K = test_age_K, 
                                      lambda1 = 1/mean_onset_time, lambda2 = 1/mean_sojourn_time))
  
  filename_tbl <- str_glue("tbl_summary_{datastamp}.csv")
  write.csv(sset, here("tables", filename_tbl))
  
}

out_sens_all <- function(preonset_rate,
                         mean_sojourn_time,
                         start_age,
                         test_age,
                         follow_up_year,
                         onset_sensitivity,
                         clinical_sensitivity,
                         dist, 
                         interval,
                         specificity,
                         saveit){
  
  dset <- expand_grid(test_age,
                     preonset_rate, 
                     MST = mean_sojourn_time,
                     start_age,
                     follow_up_year,
                     onset_sensitivity,
                     clinical_sensitivity,
                     dist, 
                     interval,
                     Specificity = specificity)
  dset <- dset %>% mutate(follow_up_age = test_age + follow_up_year)

  source(here("src", "code_analytic_sens_prosp.R"))
  dset <- dset %>% 
    rowwise() %>%
    mutate(preclinical_sens = prospective_sens_analyt(start_age = start_age,
                                                      test_age = test_age,
                                                      preonset_rate = preonset_rate,
                                                      mean_sojourn_time = MST,
                                                      onset_sensitivity = onset_sensitivity,
                                                      clinical_sensitivity = clinical_sensitivity,
                                                      dist = dist, 
                                                      interval = interval,
                                                      indolent_rate = 0,
                                                      confirmation_test_rate=1,
                                                      confirmation_test_sensitivity=1),
    ) %>%
    ungroup()
  
  source(here("src", "code_analytic_sens_retro.R"))
  dset <- dset %>% 
    rowwise() %>%
    mutate(retros_sens = analytic_retrospective_sens(t = test_age, 
                                                     t0 = start_age, 
                                                     t_check = follow_up_age,
                                                     preonset_rate = preonset_rate,
                                                     onset_sensitivity = onset_sensitivity,
                                                     clinical_sensitivity = clinical_sensitivity,
                                                     mean_sojourn_time = MST, 
                                                     dist = dist, 
                                                     interval = interval,
                                                     specificity = Specificity)) %>%
    ungroup()
  # source(here("src", "code_simulation_sens_retro.R"))
  # dset <- dset %>% 
  #   rowwise() %>%
  #   mutate(retros_sens_sim = retros_sens_sim(seed = 123, 
  #                                            nreps = 10000,
  #                                            t_check = follow_up_age,
  #                                            start_age = start_age,
  #                                            test_age = test_age,
  #                                            preonset_rate = preonset_rate,
  #                                            mean_sojourn_time = MST,
  #                                            start.dist = c(1, 0, 0),
  #                                            onset_sensitivity = onset_sensitivity,
  #                                            clinical_sensitivity = clinical_sensitivity,
  #                                            specificity = Specificity,
  #                                            pre.clinical.cancer.state = 2,
  #                                            clinical.cancer.state = 3)) %>%
  #   ungroup()
  
  source(here("src", "code_analytic_sens_em.R"))
  dset <- dset %>% 
    rowwise() %>%
    mutate(empirical_sens = analytic_empirical_sens(t = test_age, 
                                                    t0 = start_age, 
                                                    t_check = follow_up_age,
                                                    preonset_rate = preonset_rate,
                                                    onset_sensitivity = onset_sensitivity,
                                                    clinical_sensitivity = clinical_sensitivity,
                                                    mean_sojourn_time = MST,
                                                    dist = dist, 
                                                    interval = interval)) %>%
    ungroup()
  
  dset2 <- dset %>% pivot_longer(cols = c(contains("sens"), -matches("sensitivity")),#c(preclinical_sens, retros_sens, empirical_sens),
                        names_to = "Type",
                        values_to = "Sensitivity")
  dset2 <- dset2 %>% mutate(Type = ifelse(Type == "preclinical_sens", "Preclinical", 
                                          ifelse(Type == "retros_sens", "Retrospective", "Empirical")),
                            Type = factor(Type, levels = c("Preclinical", 
                                                           "Retrospective",
                                                           "Empirical")))

  if(saveit){
    filename_tbl <- str_glue("tbl_comparison_{datastamp}.csv")
    write.csv(dset, here("tables", filename_tbl))
  }
  
  return(dset2)
}

output_figures <- function(dset){
  
  gg1 <- dset2 %>% filter(Type %in% c("Preclinical", "Retrospective"),
                          test_age == 70) %>%
    ggplot(aes(x = follow_up_year, y = Sensitivity, color = Type)) +
    geom_line(size = 1.2) +
    facet_grid(MST ~ Specificity, labeller = label_both) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(limits = c(2, 20), breaks = seq(2, 20, by = 2)) +
    # scale_color_brewer(palette = "Set1") +
    labs(
      title = "Sensitivity Over Follow-up Years",
      subtitle = paste0("Start age=", start_age, ", onset rate=", preonset_rate, ", sample/screen age=", test_age[2]),
      x = "Follow-up interval (years) after blood sampling",
      y = "Sensitivity",
      color = "Type"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(color = NA, fill = NA),
      strip.text = element_text(size = 10, face = "bold"),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      strip.text.y.right = element_text(angle = 0),
      panel.spacing = unit(1, "lines"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    )
  print(gg1)
  
  gg2 <- dset2 %>% filter(Type %in% c("Preclinical", 
                                      "Empirical"),
                          Specificity == 1) %>%
    rename(`Test age` = test_age) %>%
    ggplot(aes(x = follow_up_year, y = Sensitivity, color = Type)) +
    geom_line(size = 1.2) +
    facet_grid(MST ~ `Test age`, labeller = label_both) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(limits = c(2, 20), breaks = seq(2, 20, by = 2)) +
    labs(
      title = "Sensitivity Over Follow-up Years",
      subtitle = paste0("Start age=", start_age, ", onset rate=", preonset_rate),
      x = "Follow-up interval (years) after blood sampling",
      y = "Sensitivity",
      color = "Type"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(color = NA, fill = NA),
      strip.text = element_text(size = 10, face = "bold"),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      strip.text.y.right = element_text(angle = 0),
      panel.spacing = unit(1, "lines"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    )
  
  print(gg2)
  
  if(saveit){
    filename1 <- str_glue("compare_preclinical_vs_retrospective_{datastamp}")
    ggsave(plot=print(gg1), here("figure", paste0("fig_", filename1, ".png")), width=8, height=6)
    filename2 <- str_glue("compare_preclinical_vs_empirical_{datastamp}")
    ggsave(plot=print(gg2), here("figure", paste0("fig_", filename2, ".png")), width=8, height=6)
  }
  }


control <- function(preonset_rate = 1/25,
                    mean_sojourn_time = c(2, 5, 10),
                    start_age = 40,
                    test_age = c(55, 70),
                    follow_up_year = c(1, 2, seq(4, 20, 2)),
                    onset_sensitivity = 0,
                    clinical_sensitivity = 0.8,
                    dist = "exponential", 
                    interval = c(40, 100), # for uniform distribution
                    specificity = c(0.9, 1),
                    saveit = FALSE){
  
  dset <- out_sens_all(preonset_rate,
               mean_sojourn_time,
               start_age,
               test_age,
               follow_up_year,
               onset_sensitivity,
               clinical_sensitivity,
               dist, 
               interval,
               specificity,
               saveit)
  
  summary_table(mean_onset_time = 25, 
                mean_sojourn_time = c(1, 2, 5, 10),
                threshold = 100,
                start_age = 40,
                test_age = c(55, 65, 75))
}

control(saveit = TRUE)


