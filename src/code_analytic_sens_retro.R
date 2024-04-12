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
                       is_denominator = FALSE){
  if(is_denominator){
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

f_g <- function(t_p, t_c, t0, preonset_rate, mean_sojourn_time){
  f(t_p, t0, preonset_rate) * g(t_c, t_p, mean_sojourn_time)
}

tp_integral_f_g <- function(t_c, t0, 
                            preonset_rate, mean_sojourn_time){
  integrate(Vectorize(f_g),
            lower = t0, upper = t_c,
            t_c = t_c, 
            t0 = t0, 
            preonset_rate = preonset_rate, 
            mean_sojourn_time = mean_sojourn_time)$value
}


analytic_retrospective_sens <- function(t0, t, t2,
                                      preonset_rate, mean_sojourn_time, 
                                      onset_sensitivity, clinical_sensitivity, 
                                      specificity = 0.9,
                                      method = "Estimate"){
  if(method == "Estimate"){
    true_pos <- integrate(Vectorize(integral_h_g),
                          lower = t0, upper = t,
                          t0 = t0, 
                          t = t, 
                          t2 = t2, 
                          preonset_rate = preonset_rate, 
                          mean_sojourn_time = mean_sojourn_time, 
                          onset_sensitivity = onset_sensitivity, 
                          clinical_sensitivity = clinical_sensitivity)$value
    
    false_pos <- (1-specificity) * integrate(Vectorize(integral_g),
                                             lower = t, upper = t2,
                                             t0 = t0, 
                                             t = t, 
                                             t2 = t2, 
                                             preonset_rate = preonset_rate, 
                                             mean_sojourn_time = mean_sojourn_time)$value
    
    numerator <- true_pos + false_pos
    
    denominator <- integrate(Vectorize(tp_integral_f_g),
                             lower = t, upper = t2,
                             t0 = t0, 
                             preonset_rate = preonset_rate, 
                             mean_sojourn_time = mean_sojourn_time)$value
  }
  if(method == "Truth"){
    numerator <- integrate(Vectorize(integral_h_g),
                           lower = t0, upper = t,
                           t0 = t0, 
                           t = t, 
                           t2 = t2, 
                           preonset_rate = preonset_rate, 
                           mean_sojourn_time = mean_sojourn_time, 
                           onset_sensitivity = onset_sensitivity, 
                           clinical_sensitivity = clinical_sensitivity)$value
    denominator <- integrate(Vectorize(integral_g),
                             lower = t0, upper = t,
                             t0 = t0, 
                             t = t, 
                             t2 = t2, 
                             preonset_rate = preonset_rate, 
                             mean_sojourn_time = mean_sojourn_time,
                             is_denominator = TRUE)$value
  }
  
  return(numerator/denominator)
}

