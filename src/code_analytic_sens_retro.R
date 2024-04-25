library(tidyverse)

##################################################
# Mathematical formula of retrospective sensitivity
# Author: Yibai Zhao
## t0: initial time
## t_p: pre-clinical onset time
## t_c: clinical onset time
## t2: clinical check time
## f(u): onset time ~ exp(rate)
## g(s): sojourn time ~ exp(1/mst), mst=mean sojourn time
##################################################
f <- function(u, preonset_rate){
  preonset_rate*exp(-preonset_rate*u)
}

g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}


# sensitivity function
h <- function(t, s, u, t0, onset_sensitivity, clinical_sensitivity){
  tp <- t0+u # preclinical onset time
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-tp)

  return(sensitivity)
}


# inner integral
h_g <- function(s, t, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  return(h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)*g(s, mean_sojourn_time))
}

f_integral_h_g <- function(u, t, t0, t_check, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate)*
    integrate(Vectorize(h_g),
              lower=t-tp, upper=t_check-tp,
              u=u,
              t=t,
              t0=t0,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time)$value
}

f_integral_g <- function(u, t, t0, t_check,
                       preonset_rate, mean_sojourn_time,
                       is_TRUE_positive = TRUE){
  tp <- t0+u
  if(is_TRUE_positive){
    f(u, preonset_rate)*
      integrate(Vectorize(g),
                lower = t-tp, upper = t_check-tp,
                mean_sojourn_time = mean_sojourn_time)$value
  }else{
    f(u, preonset_rate)*
      integrate(Vectorize(g),
                lower = 0, upper = t_check-tp,
                mean_sojourn_time = mean_sojourn_time)$value
  }
}

integral_f_h_g <- function(t, t0, t_check, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  integrate(Vectorize(f_integral_h_g),
            lower = 0, upper = t-t0,
            t = t,
            t0 = t0, 
            t_check = t_check, 
            preonset_rate = preonset_rate, 
            onset_sensitivity = onset_sensitivity, 
            clinical_sensitivity = clinical_sensitivity,
            mean_sojourn_time = mean_sojourn_time)$value
}

integral_f_g <- function(t, t0, t_check,
                         preonset_rate, mean_sojourn_time,
                         is_TRUE_positive = TRUE){
  if(is_TRUE_positive){
    integrate(Vectorize(f_integral_g),
              lower = 0, upper = t-t0,
              t = t,
              t0 = t0, 
              t_check = t_check, 
              preonset_rate = preonset_rate, 
              mean_sojourn_time = mean_sojourn_time,
              is_TRUE_positive = TRUE)$value
  }else{
    integrate(Vectorize(f_integral_g),
              lower = t-t0, upper = t_check-t0,
              t = t,
              t0 = t0, 
              t_check = t_check, 
              preonset_rate = preonset_rate, 
              mean_sojourn_time = mean_sojourn_time,
              is_TRUE_positive = FALSE)$value
    
  }
}


analytic_retrospective_sens <- function(t, t0, t_check,
                                        preonset_rate,
                                        onset_sensitivity,
                                        clinical_sensitivity,
                                        mean_sojourn_time, 
                                        specificity = 1){
  
  numerator_TRUE_positive <- integral_f_h_g(t, t0, t_check,
                                          preonset_rate, 
                                          onset_sensitivity, 
                                          clinical_sensitivity,
                                          mean_sojourn_time)
  numerator_FALSE_positive <- (1 - specificity) *
    integral_f_g(t, t0, t_check,
                 preonset_rate, mean_sojourn_time,
                 is_TRUE_positive = FALSE)
  
  denominator_TRUE_positive <- integral_f_g(t, t0, t_check,
                                            preonset_rate, mean_sojourn_time,
                                            is_TRUE_positive = TRUE)
  denominator_FALSE_positive <- integral_f_g(t, t0, t_check,
                                            preonset_rate, mean_sojourn_time,
                                            is_TRUE_positive = FALSE)
  
  return((numerator_TRUE_positive + numerator_FALSE_positive) / (denominator_TRUE_positive + denominator_FALSE_positive))
}

