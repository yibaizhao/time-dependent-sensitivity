library(tidyverse)
library(ggplot2)

# datestamp <- "2023-11-30"
datestamp <- "2023-12-05"

set.seed(1234)

##################################################
# Prospective sensitivity evaluated at specified times
# sojourn_time: sojourn time in years
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
# TODO: Determine annotation positions programmatically
##################################################
test_sensitivity <- function(sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             is_indolent=FALSE,
                             time=seq(0, 20),
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
    tset <- tset %>%
      mutate(sensitivity=ifelse(sensitivity<=clinical_sensitivity, sensitivity, NA),
             sensitivity=ifelse(time<0|time>sojourn_time, NA, sensitivity))
  }else{
    tset <- tibble(time=time, 
                   sensitivity=onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/sojourn_time*time)))
  }
  return(tset)
}

##################################################
# Simulate natural histories and determine sensitivity
# of a test at the specified age that would be estimated
# in a retrospective setting among persons with clinical
# cancer
# N: number of simulated individuals
# start_age: initial age at start of simulation
# preonset_rate: exponential rate of cancer onset
# mean_sojourn_time: exponential mean sojourn time
# indolent_rate: proportion of non-progressive cancers
# test_age: age at one-time screening test
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
##################################################
retrospective_test_sensitivity <- function(N,
                                         preonset_rate,
                                         mean_sojourn_time,
                                         sample_window,
                                         start_age,
                                         retro_age,
                                         onset_sensitivity,
                                         clinical_sensitivity,
                                         specificity){
  sample_age = retro_age - sample_window
  # simulate ages at preclinical onset and clinical diagnosis
  dset <- tibble(id = 1:N,
                 onset_age = start_age + rexp(N, rate = preonset_rate),
                 sojourn_time = rexp(N, rate = 1/mean_sojourn_time),
                 clinical_age = onset_age + sojourn_time)
  # find cases that clinically diagnosed at test age
  cset <- dset %>% filter(clinical_age <= retro_age & clinical_age > sample_age)
  # estimate screen sensitivity at sample time
  cset <- cset %>%
    mutate(screen_sens = test_sensitivity(sojourn_time = sojourn_time,
                                       onset_sensitivity,
                                       clinical_sensitivity,
                                       is_indolent = FALSE,
                                       time = sample_age-onset_age)$sensitivity)
  # if clinical diagnosed before sample time, then sensitivity is 80% (need to verify)
  cset <- cset %>% mutate(screen_sens = ifelse(clinical_age <= sample_age, clinical_sensitivity, screen_sens))
  # if preclinical onset after sample time, then sensitivity is 0% (need to verify)
  cset <- cset %>% mutate(screen_sens = ifelse(onset_age > sample_age, 0, screen_sens))
  # estimate test results at sample time
  cset <- cset %>% 
    group_by(id) %>%
    mutate(test_results = ifelse(onset_age <= sample_age, rbinom(1, 1, screen_sens), rbinom(1, 1, 1-!!specificity))) %>%
    ungroup()
  
  retros_sens <- cset %>%
    summarise(retros_sens = mean(test_results))
  return(retros_sens$retros_sens)
}


plot_prospective_sensitivity <- function(N,
                                         preonset_rate,
                                         mean_sojourn_time_seq,
                                         sample_window_seq,
                                         start_age,
                                         retro_age,
                                         onset_sensitivity,
                                         clinical_sensitivity,
                                         specificity = 1){
  sset <- data.frame(N,
                     start_age,
                     preonset_rate,
                     retro_age)
  sset <- expand_grid(sset, 
                      sample_window = sample_window_seq,
                      mean_sojourn_time = mean_sojourn_time_seq)
  sset <- sset %>%
    group_by(mean_sojourn_time, sample_window) %>%
    mutate(retro_sens = retrospective_test_sensitivity(N = N,
                                                       preonset_rate = preonset_rate,
                                                       mean_sojourn_time = mean_sojourn_time,
                                                       sample_window = sample_window,
                                                       start_age = start_age,
                                                       retro_age = retro_age,
                                                       onset_sensitivity = !!onset_sensitivity,
                                                       clinical_sensitivity = !!clinical_sensitivity,
                                                       specificity = !!specificity),
           retro_sens_analytic = analytic_prospective_sens(t0 = start_age, 
                                                           t = retro_age-sample_window, 
                                                           t2 = retro_age, 
                                                           preonset_rate = preonset_rate,
                                                           mean_sojourn_time = mean_sojourn_time, 
                                                           onset_sensitivity = !!onset_sensitivity, 
                                                           clinical_sensitivity = !!clinical_sensitivity,
                                                           specificity = !!specificity)) %>%
    ungroup()
  
  sset <- sset %>% pivot_longer(cols = c("retro_sens", "retro_sens_analytic"),
                                values_to = "Retrospective sensitivity",
                                names_to = "Method")
  sset %>%
    ggplot(aes(x = retro_age - sample_window, y = `Retrospective sensitivity`, color = Method)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = retro_age-sample_window_seq) +    
    facet_grid(.~mean_sojourn_time) +
    xlab("Sample age = Retro age (60) - sample window") +
    ylab("Retrospective sensitivity") +
    ylim(0, 1)
}

##################################################
# Mathematical formula of retrospective sensitivity
## t0: initial time
## t_p: pre-clinical onset time
## t_c: clinical onset time
## t2: clinical check time
## f(u): onset time ~ exp(rate)
## g(s): sojourn time ~ exp(1/mst), mst=mean sojourn time
##################################################
f <- function(t_p, t0, preonset_rate){
  u <- t_p-t0
  preonset_rate*exp(-preonset_rate*u)
}

g <- function(t_c, t_p, mean_sojourn_time){
  s <- t_c-t_p
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

# true screening sensitivity function
h <- function(t_c, t_p, t, onset_sensitivity, clinical_sensitivity){
  # sojourn time
  s <- t_c-t_p
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-t_p)
  sensitivity <- ifelse(sensitivity<=clinical_sensitivity, sensitivity, NA)
  
  return(sensitivity)
}

# inner integral
h_g <- function(t_c, t_p, t, mean_sojourn_time, onset_sensitivity, clinical_sensitivity){
  h(t_c, t_p, t, onset_sensitivity, clinical_sensitivity) * g(t_c, t_p, mean_sojourn_time)
}

integral_h_g <- function(t_p, t0, t, t2, 
                         preonset_rate, mean_sojourn_time, 
                         onset_sensitivity, clinical_sensitivity){
  f(t_p, t0, preonset_rate)*
    integrate(Vectorize(h_g),
              lower = t, upper = t2,
              t_p=t_p, 
              t = t, 
              mean_sojourn_time = mean_sojourn_time, 
              onset_sensitivity = onset_sensitivity, 
              clinical_sensitivity = clinical_sensitivity)$value
}

integral_g <- function(t_p, t0, t, t2, 
                      preonset_rate, mean_sojourn_time,
                      denominator = FALSE){
  if(denominator){
    f(t_p, t0, preonset_rate)*
      integrate(Vectorize(g),
                lower = t, upper = t2,
                t_p = t_p, 
                mean_sojourn_time = mean_sojourn_time)$value
  }else{
    f(t_p, t0, preonset_rate)*
      integrate(Vectorize(g),
                lower = t_p, upper = t2,
                t_p = t_p, 
                mean_sojourn_time = mean_sojourn_time)$value
  }
}

analytic_prospective_sens <- function(t0, t, t2,
                                      preonset_rate, mean_sojourn_time, 
                                      onset_sensitivity, clinical_sensitivity, 
                                      specificity = 0.9){
  true_pos <- integrate(Vectorize(integral_h_g),
                        lower = t0, upper = t,
                        t0 = t0, 
                        t = t, 
                        t2 = t2, 
                        preonset_rate = preonset_rate, 
                        mean_sojourn_time = mean_sojourn_time, 
                        onset_sensitivity = onset_sensitivity, 
                        clinical_sensitivity = clinical_sensitivity)$value
  
  false_pos <- specificity * integrate(Vectorize(integral_g),
                                       lower = t, upper = t2,
                                       t0 = t0, 
                                       t = t, 
                                       t2 = t2, 
                                       preonset_rate = preonset_rate, 
                                       mean_sojourn_time = mean_sojourn_time)$value
  
  numerator <- true_pos + false_pos
  
  denominator <- integrate(Vectorize(integral_g),
                            lower = t0, upper = t2,
                            t0 = t0, 
                            t = t, 
                            t2 = t2, 
                            preonset_rate = preonset_rate, 
                            mean_sojourn_time = mean_sojourn_time,
                            denominator = TRUE)$value
  
  return(numerator/denominator)
}

control <- function(N=1000,
                    preonset_rate=0.1,
                    sample_window_seq = seq(2, (60-40-2), by = 2),
                    mean_sojourn_time_seq = c(2, 4, 6, 10),
                    start_age=40,
                    retro_age=60,
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    specificity=0.9
){
  plot_prospective_sensitivity(N,
                               preonset_rate,
                               mean_sojourn_time_seq,
                               sample_window_seq,
                               start_age,
                               retro_age,
                               onset_sensitivity,
                               clinical_sensitivity,
                               specificity)
}

control()

