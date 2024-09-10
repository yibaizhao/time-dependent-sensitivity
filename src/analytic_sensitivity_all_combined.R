library(here)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(scales)
library(foreach)
library(doParallel)
library(magrittr)
library(readr)

# datastamp <- "2024-04-15"
# datastamp <- "2024-04-17"
# datastamp <- "2024-04-24"
# datastamp <- "2024-04-25"
# datastamp <- "2024-05-01"
# datastamp <- "2024-05-07"
# datastamp <- "2024-05-01"
# datastamp <- "2024-07-15"
# datastamp <- "2024-07-16"
# datastamp <- "2024-07-17"
# datastamp <- "2024-07-29"
# datastamp <- "2024-08-13"
datastamp <- "2024-09-10"

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

# Prob(tp<=test_age_range<tc)
CDF_tp_K_tc_all_age <- function(t0, test_age_range_lb, test_age_range_ub, lambda1, lambda2){
  integral_all_age <- integrate(Vectorize(CDF_tp_K_tc),
            lower = test_age_range_lb, upper = test_age_range_ub,
            t0 = t0, 
            lambda1 = lambda1, 
            lambda2 = lambda2)$value
  return(integral_all_age / (test_age_range_ub - test_age_range_lb))
}

# calculate t_p<=C and t_c<=C
summary_table <- function(mean_onset_time, 
                          mean_sojourn_time,
                          threshold,
                          start_age,
                          test_age_range){
  lambda1 <- 1/mean_onset_time # parameter for onset time ~ exp(lambda1)
  lambda2 <- 1/mean_sojourn_time # parameter for sojourn time ~ exp(lambda2)

  sset <- expand.grid(start_age = start_age,
                      test_age_range_lb = test_age_range[1],
                      test_age_range_ub = test_age_range[2],
                      mean_onset_time = mean_onset_time, 
                      mean_sojourn_time = mean_sojourn_time,
                      threshold_C = threshold)
  
  sset <- sset %>% 
    mutate(prob_tp_C = CDF_tp(t0 = start_age, C = threshold_C, 
                              lambda1 = 1/mean_onset_time),
           prob_tc_C = CDF_tp_tc(t0 = start_age, C = threshold_C, 
                                 lambda1 = 1/mean_onset_time, lambda2 = 1/mean_sojourn_time),
           CDF_tp_K_tc_all_age = CDF_tp_K_tc_all_age(t0 = start_age, 
                                             test_age_range_lb = test_age_range_lb, 
                                             test_age_range_ub = test_age_range_ub, 
                                             lambda1 = 1/mean_onset_time, 
                                             lambda2 = 1/mean_sojourn_time))
  
  filename_tbl <- str_glue("tbl_summary_{datastamp}.csv")
  write.csv(sset, here("tables", filename_tbl))
  
}

out_sens_all <- function(preonset_rate,
                         mean_sojourn_time,
                         shift,
                         start_age,
                         test_age_range, # 50-74
                         follow_up_year,
                         onset_sensitivity,
                         clinical_sensitivity,
                         dist, 
                         interval,
                         specificity,
                         saveit){
  
  dset <- expand_grid(preonset_rate, 
                     MST = mean_sojourn_time,
                     start_age,
                     test_age_range_lb = test_age_range[1],
                     test_age_range_ub = test_age_range[2],
                     follow_up_year,
                     onset_sensitivity,
                     clinical_sensitivity,
                     dist, 
                     interval,
                     Specificity = specificity)
  # dset <- dset %>% mutate(follow_up_age = test_age + follow_up_year)
  dset <- dset %>% mutate(row_id = 1:nrow(dset))

  # Initialize parallel backend
  registerDoParallel(cores = detectCores())
  
  # Define the function to process each subset of dset
  process_subset <- function(dset_i) {
    # Apply transformations and calculate sensitivities
    source(here("src", "code_analytic_sens_prosp.R"))
    dset_i <- dset_i %>%
      rowwise() %>%
      mutate(
        preclinical_sens = prospective_sens_all_test_age_analyt(
          test_age_range = c(test_age_range_lb, test_age_range_ub),
          start_age = start_age,
          shift = shift,
          preonset_rate = preonset_rate,
          mean_sojourn_time = MST,
          onset_sensitivity = onset_sensitivity,
          clinical_sensitivity = clinical_sensitivity,
          dist = dist, 
          interval = interval,
          indolent_rate = 0,
          confirmation_test_rate = 1,
          confirmation_test_sensitivity = 1
        ))
    source(here("src", "code_analytic_sens_retro.R"))
    dset_i <- dset_i %>%
      rowwise() %>%
      mutate(
        retros_sens = analytic_retrospective_sens_all_test_age(
          test_age_range = c(test_age_range_lb, test_age_range_ub),
          t0 = start_age, 
          shift = shift, 
          follow_up_year = follow_up_year,
          preonset_rate = preonset_rate,
          onset_sensitivity = onset_sensitivity,
          clinical_sensitivity = clinical_sensitivity,
          mean_sojourn_time = MST, 
          dist = dist, 
          interval = interval,
          specificity = Specificity
        ))
    source(here("src", "code_analytic_sens_em.R"))
    dset_i <- dset_i %>%
      rowwise() %>%
      mutate(
        empirical_sens = analytic_empirical_sens_all_test_age(
          test_age_range = c(test_age_range_lb, test_age_range_ub),
          t0 = start_age, 
          shift = shift, 
          follow_up_year = follow_up_year,
          preonset_rate = preonset_rate,
          onset_sensitivity = onset_sensitivity,
          clinical_sensitivity = clinical_sensitivity,
          mean_sojourn_time = MST,
          dist = dist, 
          interval = interval
        )
      ) %>%
      ungroup()
    
    return(dset_i)
  }
  
  # Use parallel execution with %dopar%
  dset <- foreach(i = unique(dset$row_id), .combine = rbind, .packages = c("dplyr", "here")) %dopar% {
    dset_i <- dset %>% filter(row_id == i)
    processed_subset <- process_subset(dset_i)
    return(processed_subset)
  }
  
  dset2 <- dset %>% pivot_longer(cols = c(contains("sens"), -matches("sensitivity")),#c(preclinical_sens, retros_sens, empirical_sens),
                        names_to = "Type",
                        values_to = "Sensitivity")
  dset2 <- dset2 %>% mutate(Type = ifelse(Type == "preclinical_sens", "Preclinical", 
                                          ifelse(Type == "retros_sens", "Retrospective", "Empirical")),
                            Type = factor(Type, levels = c("Preclinical", 
                                                           "Retrospective",
                                                           "Empirical")))

  if(saveit){
    filename_tbl <- str_glue("tbl_all_sens_{datastamp}.csv")
    write.csv(dset2, here("tables", filename_tbl))
  }
  
  return(dset2)
}

output_figures <- function(dset){
  dset <- dset %>% filter(clinical_sensitivity == 0.8)
  
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))

  # Figure 1: Phase I, mean preclinical sensitivity across MST
  dset_preclinical <- dset %>%
    select(-c(row_id, follow_up_year, interval, Specificity)) %>%
    filter(Type %in% c("Preclinical")) %>%
    unique() %>%
    mutate(degradation = clinical_sensitivity - Sensitivity,
           Sensitivity = percent(Sensitivity),
           degradation = percent(degradation))
  # Define colors for each level of mean_sojourn_time
  colors <- c("1" = "grey10", "2" = "grey30", "4" = "grey50", "6" = "grey70")
  
  gg1 <- dset_preclinical %>%
    ggplot(aes(x = factor(MST), y = Sensitivity, fill = factor(MST))) +
    geom_hline(aes(yintercept = clinical_sensitivity), linetype = 'dashed') +
    geom_bar(stat = "identity", width = 1/2, position = position_dodge(width = 0.4)) +
    geom_linerange(aes(ymin = Sensitivity, ymax = clinical_sensitivity, color = factor(MST)),
                   alpha = 0.3, position = position_dodge2(width = 1)) +
    geom_text(aes(label = sprintf('-%2.0f%%', 100 * degradation), 
                  y = Sensitivity + degradation / 2, color = factor(MST)),
              size = 3, position = position_dodge2(width = 1), show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    # facet_grid(. ~ clinical_sensitivity, labeller = label_both) +
    labs(
      y = "Mean preclinical sensitivity",
      x = "Mean sojourn time"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = label_percent(accuracy = 1),
      expand = c(0, 0)
    ) +
    theme(legend.position = "none")  
  print(gg1)
    
  dset_preclinical_retrospective <- dset %>%
    select(-c(row_id, interval)) %>%
    filter(Type %in% c("Preclinical", "Retrospective"),
           # MST %in% c(2, 5, 10),
           follow_up_year <= 10) %>%
    mutate(Type = factor(Type, levels = c("Preclinical", "Retrospective"))) %>%
    unique()
  gg2 <- dset_preclinical_retrospective %>%
    ggplot(aes(x = follow_up_year, y = Sensitivity, color = Type, linetype = factor(Specificity))) +
    geom_line(size = 1.2) +
    facet_grid(MST ~., labeller = label_both) +
    scale_y_continuous(limits=c(0, 1),
                       breaks=seq(0, 1, by=0.2),
                       labels=label_percent(accuracy=1),
                       expand=c(0, 0)) +
    scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, by = 1)) +
    # scale_color_brewer(palette = "Set1") +
    labs(
      # title = "Sensitivity Over Follow-up Years",
      # subtitle = paste0("Start age=", start_age, ", onset rate=", preonset_rate,
      #                   ", clinical sensitivity=", unique(dset_preclinical_retrospective$clinical_sensitivity)),
      x = "Follow-up interval (years) after blood sampling",
      y = "Sensitivity",
      color = "Type",
      linetype = "Specificity"
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
  
  # Figure 3: Phase III, preclincal VS empirical sensitivity
  dset_empirical_Jane <- data.frame(MST = c(1, 2, 3, 4, 6, 1, 2, 3, 4, 6),
                                    follow_up_year = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                                    Sensitivity = c(0.35, 0.5, 0.61, 0.7, 0.77, 0.2, 0.33, 0.43, 0.5, 0.6),
                                    Type = rep("Empirical", 10))
  dset_preclinical_empirical <- dset %>% 
    select(-c(row_id, interval, Specificity)) %>%
    filter(Type %in% c("Preclinical", "Empirical"),
           follow_up_year %in% c(1, 2),
           clinical_sensitivity == 0.8
           # MST %in% c(1, 2, 3, 5)
           ) %>%
    mutate(Type = factor(Type, levels = c("Preclinical", "Empirical"))) %>%
    unique()
  gg3 <- dset_preclinical_empirical %>%
    ggplot(aes(x = MST, y = Sensitivity, color = Type, shape = Type, linetype = factor(follow_up_year))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    # geom_point(data = dset_empirical_Jane, aes(x = MST, y = Sensitivity, color = Type, shape = Type), color = "yellow") +
    # facet_grid(.~ follow_up_year, labeller = label_both) +
    scale_y_continuous(limits=c(0, 1),
                       breaks=seq(0, 1, by=0.2),
                       labels=label_percent(accuracy=1),
                       expand=c(0, 0)) +
    scale_x_continuous(breaks = unique(dset_preclinical_empirical$MST)) +
    labs(
      # title = "Sensitivity over Mean Sojourn Time",
      # subtitle = paste0("Start age=", start_age, ", onset rate=", preonset_rate,
      #                   ", clinical sensitivity=", unique(dset_preclinical_empirical$clinical_sensitivity)),
      x = "Mean Sojourn Time",
      y = "Sensitivity",
      color = "Type",
      linetype = "Followup Year"
    ) +
    theme_classic(base_size = 12)+
    guides(
      linetype = guide_legend(order = 2),
      color = guide_legend(order = 1),
      shape = guide_legend(order = 1)
    ) +
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
  
  print(gg3)
  
  if(saveit){
    filename1 <- str_glue("compare_preclinical_across_MST_{datastamp}")
    ggsave(plot=print(gg1), here("figure", paste0("fig_", filename1, ".pdf")), width=6, height=6)
    filename2 <- str_glue("compare_preclinical_vs_retrospective_{datastamp}")
    ggsave(plot=print(gg2), here("figure", paste0("fig_", filename2, ".pdf")), width=8, height=5)
    filename3 <- str_glue("compare_preclinical_vs_empirical_{datastamp}")
    ggsave(plot=print(gg3), here("figure", paste0("fig_", filename3, ".pdf")), width=8, height=5)
  }
  }


control <- function(preonset_rate = 1/25,
                    mean_sojourn_time = c(1, 2, 4, 6),
                    start_age = 40,
                    test_age_range = c(50, 74),
                    shift = 1,
                    follow_up_year = c(1:10),
                    onset_sensitivity = 0,
                    clinical_sensitivity = c(0.8, 0.95),
                    dist = "exponential", 
                    interval = NA, # for uniform distribution
                    specificity = c(0.9, 1),
                    saveit = TRUE){
  
  out_sens_all(preonset_rate,
               mean_sojourn_time,
               start_age,
               test_age_range,
               shift,
               follow_up_year,
               onset_sensitivity,
               clinical_sensitivity,
               dist,
               interval,
               specificity,
               saveit = saveit)
  # dset <- read_csv(here("tables", "tbl_all_sens_2024-07-17.csv"))[,-1]
  # dset <- read_csv(here("tables", "tbl_all_sens_2024-08-13.csv"))[,-1]
  # summary_table(mean_onset_time = 25, 
  #               mean_sojourn_time = c(1, 2, 3, 5, 10),
  #               threshold = 100,
  #               start_age = 40,
  #               test_age_range = c(50, 74))
}

control(saveit = TRUE)


