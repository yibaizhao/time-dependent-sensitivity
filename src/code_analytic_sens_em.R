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
f <- function(u, preonset_rate=NULL, dist="exponential", interval=NA){
  if(dist=="exponential"){
    return(preonset_rate*exp(-preonset_rate*u))
  }
  if(dist=="uniform"){
    return(1/(interval[2]-interval[1]))
  }
}

# sojourn time distribution
g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

h_g <- function(t, u, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  h(t, u, s, t0, onset_sensitivity, clinical_sensitivity)*g(s, mean_sojourn_time)
}

f_integral_h_g <- function(u, t, t0, shift, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time,
                           dist, interval){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate, dist, interval)*
    integrate(Vectorize(h_g),
              lower = max(t-tp, shift), upper = Inf,
              t = t,
              u = u,
              t0 = t0,
              onset_sensitivity = onset_sensitivity,
              clinical_sensitivity = clinical_sensitivity,
              mean_sojourn_time = mean_sojourn_time)$value
}

screen_detected <- function(t, t0, shift,
                            preonset_rate,
                            onset_sensitivity,
                            clinical_sensitivity,
                            mean_sojourn_time,
                            dist, 
                            interval){
  integrate(Vectorize(f_integral_h_g),
            lower = 0, upper = t-t0,
            t = t,
            t0 = t0, 
            shift = shift,
            preonset_rate = preonset_rate,
            onset_sensitivity = onset_sensitivity,
            clinical_sensitivity = clinical_sensitivity,
            mean_sojourn_time = mean_sojourn_time,
            dist = dist, 
            interval = interval)$value
}

oppo_h_g <- function(t, u, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  (1 - h(t, u, s, t0, onset_sensitivity, clinical_sensitivity)) * g(s, mean_sojourn_time)
}

f_integral_oppo_h_g <- function(u, t, t0, t_check, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time, 
                                dist, interval){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate, dist, interval)*
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
                                           mean_sojourn_time,
                                           dist, interval){
  integrate(Vectorize(f_integral_oppo_h_g),
            lower = 0, upper = t-t0,
            t = t,
            t0 = t0,
            t_check = t_check, 
            preonset_rate = preonset_rate,
            onset_sensitivity = onset_sensitivity,
            clinical_sensitivity = clinical_sensitivity,
            mean_sojourn_time = mean_sojourn_time,
            dist = dist, 
            interval = interval)$value
}

f_integral_g <- function(u, t0, t_check, shift, preonset_rate, mean_sojourn_time, 
                         dist, interval){
  tp <- t0+u # preclinical onset time
  f(u, preonset_rate, dist, interval)*
    integrate(Vectorize(g),
              lower = shift, upper = t_check - tp,
              mean_sojourn_time = mean_sojourn_time)$value
}

new_case <- function(t, t0, t_check, shift,
                     preonset_rate,
                     mean_sojourn_time,
                     dist, interval){
  integrate(Vectorize(f_integral_g),
            lower = t - t0, upper = t_check - t0,
            t0 = t0,
            t_check = t_check, 
            shift = shift, 
            preonset_rate = preonset_rate,
            mean_sojourn_time = mean_sojourn_time,
            dist = dist, 
            interval = interval)$value
}

analytic_empirical_sens <- function(t, t0, shift, follow_up_year,
                                    preonset_rate,
                                    onset_sensitivity,
                                    clinical_sensitivity,
                                    mean_sojourn_time,
                                    dist, 
                                    interval){
  t_check <- t + follow_up_year
  screen_detected <- screen_detected(t, t0, shift, 
                               preonset_rate,
                               onset_sensitivity,
                               clinical_sensitivity,
                               mean_sojourn_time,
                               dist, 
                               interval)
  interval_cancer_false_negative <- interval_cancer_false_negative(t, t0, t_check,
                                                                preonset_rate,
                                                                onset_sensitivity,
                                                                clinical_sensitivity,
                                                                mean_sojourn_time,
                                                                dist, interval)
  interval_cancer_new_case <- new_case(t, t0, t_check, shift,
                                       preonset_rate,
                                       mean_sojourn_time,
                                       dist, interval)
  
  return(screen_detected / (screen_detected + interval_cancer_false_negative + interval_cancer_new_case))
}

analytic_empirical_sens_all_test_age <- function(test_age_range,
                                                 t0, shift, follow_up_year,
                                                  preonset_rate,
                                                  onset_sensitivity,
                                                  clinical_sensitivity,
                                                  mean_sojourn_time,
                                                  dist, 
                                                  interval){
  overall_sensitivity <- integrate(Vectorize(analytic_empirical_sens),
                                   lower=test_age_range[1], upper=test_age_range[2],
                                   t0 = t0, 
                                   shift = shift, 
                                   follow_up_year = follow_up_year,
                                   preonset_rate = preonset_rate,
                                   onset_sensitivity = onset_sensitivity,
                                   clinical_sensitivity = clinical_sensitivity,
                                   mean_sojourn_time = mean_sojourn_time,
                                   dist = dist, 
                                   interval = interval)$value
  
  return(overall_sensitivity / (test_age_range[2] - test_age_range[1]))
  
}
