library(here)
library(ggplot2)
library(stringr)
library(tidyr)

# datastamp <- "2024-04-15"
datastamp <- "2024-04-17"

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

out_sens_all <- function(preonset_rate,
                         mean_sojourn_time,
                         start_age,
                         test_age,
                         follow_up_year,
                         onset_sensitivity,
                         clinical_sensitivity,
                         specificity,
                         saveit){
  
  dset <- expand_grid(test_age,
                     preonset_rate, 
                     MST = mean_sojourn_time,
                     start_age,
                     follow_up_year,
                     onset_sensitivity,
                     clinical_sensitivity,
                     Specificity = specificity)
  dset$follow_up_age <- dset$test_age+dset$follow_up_year

  source(here("src", "code_analytic_sens_prosp.R"))
  dset <- dset %>% 
    rowwise() %>%
    mutate(preclinical_sens = prospective_sens_analyt(start_age = start_age,
                                                     test_age = test_age,
                                                     preonset_rate = preonset_rate,
                                                     mean_sojourn_time = MST,
                                                     onset_sensitivity = onset_sensitivity,
                                                     clinical_sensitivity = clinical_sensitivity,
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
                                                     specificity = Specificity)) %>%
    ungroup()
  
  dset2 <- dset %>% pivot_longer(cols = c(preclinical_sens, retros_sens),
                        names_to = "Type",
                        values_to = "Sensitivity")
  dset2 <- dset2 %>% mutate(Type = ifelse(Type == "preclinical_sens", "Pre-clinical", "Retrospective"))
  # figure
  gg <- ggplot(dset2, aes(x = factor(follow_up_year), y = Sensitivity, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    facet_grid(MST ~ Specificity, 
               labeller = label_both) +
    ylim(0, 1) +
    theme(legend.position = "top") +
    labs(subtitle = paste0("Start age=", start_age, ", onset rate=", preonset_rate, ", sample/screen age=", test_age),
         x = " Years follow-up after blood sample")
  print(gg)
  
  if(saveit){
    filename <- str_glue("compare_preclinical_vs_retrospective_{datastamp}")
    ggsave(plot=print(gg), here("figure", paste0("fig_", filename, ".png")), width=8, height=6)
    write.csv(dset, here("tables", paste0("tbl_", filename, ".csv")))
  }
}


control <- function(preonset_rate = 0.1,
                    mean_sojourn_time = c(5, 10, 20),
                    start_age = 40,
                    test_age = 50,
                    follow_up_year = c(5, 10, 20),
                    onset_sensitivity = 0.3,
                    clinical_sensitivity = 0.8,
                    specificity = c(0.9, 1),
                    saveit = FALSE){
  
  out_sens_all(preonset_rate,
               mean_sojourn_time,
               start_age,
               test_age,
               follow_up_year,
               onset_sensitivity,
               clinical_sensitivity,
               specificity,
               saveit)
}

control(saveit = TRUE)




source(here("src", "code_analytic_sens_em.R"))
rate.matrix <- matrix(c(-preonset_rate, preonset_rate, 0,
                        0, -1/mean_sojourn_time, 1/mean_sojourn_time,
                        0, 0, 0),
                      byrow=T,nrow=3)
sens_em <- sapply(test_age, function(t) {
  empirical.sensitivity.simulation(seed=123, 
                                   rate.matrix = rate.matrix,
                                   start.dist = c(1,0,0),
                                   start_age = start_age,
                                   screen.times = t,
                                   post.screen.lookout = screen_interval,
                                   clinical.cancer.state = 3,
                                   pre.clinical.cancer.state = 2,
                                   onset_sensitivity = onset_sensitivity,
                                   clinical_sensitivity = clinical_sensitivity,
                                   multiple.tests = FALSE,
                                   nreps = 10000)
  })
sens_em <- t(as.matrix(sens_em)) %>% as.data.frame() %>% unnest()
sens_em$test_age <- test_age
sens_em$type <- "Empirical"

