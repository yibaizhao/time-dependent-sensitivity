##################################################
# Mathematical formula of empirical sensitivity
# Author: Yibai Zhao
##################################################
# sensitivity function
h <- function(t, u, s, t0, onset_sensitivity, clinical_sensitivity){
  tp <- t0+u # preclinical onset time
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-tp)

  return(sensitivity)
}

# onset distribution
f <- function(u, preonset_rate){
  preonset_rate*exp(-preonset_rate*u)
}

# sojourn time distribution
g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

h_g <- function(t, u, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  h(t, u, s, t0, onset_sensitivity, clinical_sensitivity)*g(s, mean_sojourn_time)
}

f_integral_h_g <- function(u, t, t0, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate)*
    integrate(Vectorize(h_g),
              lower = t-tp, upper = Inf,
              t = t,
              u = u,
              t0 = t0,
              onset_sensitivity = onset_sensitivity,
              clinical_sensitivity = clinical_sensitivity,
              mean_sojourn_time = mean_sojourn_time)$value
}

screen_detected <- function(t, t0,
                            preonset_rate,
                            onset_sensitivity,
                            clinical_sensitivity,
                            mean_sojourn_time){
  integrate(Vectorize(f_integral_h_g),
            lower = 0, upper = t-t0,
            t = t,
            t0 = t0,
            preonset_rate = preonset_rate,
            onset_sensitivity = onset_sensitivity,
            clinical_sensitivity = clinical_sensitivity,
            mean_sojourn_time = mean_sojourn_time)$value
}

oppo_h_g <- function(t, u, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  (1 - h(t, u, s, t0, onset_sensitivity, clinical_sensitivity)) * g(s, mean_sojourn_time)
}

f_integral_oppo_h_g <- function(u, t, t0, t_check, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate)*
    integrate(Vectorize(oppo_h_g),
              lower = t-tp, upper = t_check-tp,
              t = t,
              u = u,
              t0 = t0,
              onset_sensitivity = onset_sensitivity,
              clinical_sensitivity = clinical_sensitivity,
              mean_sojourn_time = mean_sojourn_time)$value
}

interval_cancer_false_negative <- function(t, t0, t_check,
                                           preonset_rate,
                                           onset_sensitivity,
                                           clinical_sensitivity,
                                           mean_sojourn_time){
  integrate(Vectorize(f_integral_oppo_h_g),
            lower = 0, upper = t-t0,
            t = t,
            t0 = t0,
            t_check = t_check, 
            preonset_rate = preonset_rate,
            onset_sensitivity = onset_sensitivity,
            clinical_sensitivity = clinical_sensitivity,
            mean_sojourn_time = mean_sojourn_time)$value
}

f_integral_g <- function(u, t0, t_check, preonset_rate, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate)*
    integrate(Vectorize(g),
              lower = 0, upper = t_check - tp,
              mean_sojourn_time = mean_sojourn_time)$value
}

new_case <- function(t, t0, t_check,
                     preonset_rate,
                     mean_sojourn_time){
  integrate(Vectorize(f_integral_g),
            lower = t - t0, upper = t_check - t0,
            t0 = t0,
            t_check = t_check, 
            preonset_rate = preonset_rate,
            mean_sojourn_time = mean_sojourn_time)$value
}

analytic_empirical_sens <- function(t, t0, t_check,
                                    preonset_rate,
                                    onset_sensitivity,
                                    clinical_sensitivity,
                                    mean_sojourn_time){
  screen_detected <- screen_detected(t, t0, 
                               preonset_rate,
                               onset_sensitivity,
                               clinical_sensitivity,
                               mean_sojourn_time)
  interval_cancer_false_negative <- interval_cancer_false_negative(t, t0, t_check,
                                                                preonset_rate,
                                                                onset_sensitivity,
                                                                clinical_sensitivity,
                                                                mean_sojourn_time)
  interval_cancer_new_case <- new_case(t, t0, t_check,
                       preonset_rate,
                       mean_sojourn_time)
  
  return(screen_detected / (screen_detected + interval_cancer_false_negative + interval_cancer_new_case))
}
