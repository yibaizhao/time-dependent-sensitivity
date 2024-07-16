library(tidyverse)

##################################################
# Mathematical formula of prospective sensitivity
# Author: Yibai Zhao
## f(u, preonset_rate): function of preclinical onset time and onset rate
##            u: time interval between initial age to preclinical onset
## g(s, mean_sojourn_time): function of sojourn time
##            s: sojourn time, time between preclinical onset and clinical onset
##            mean_sojourn_time: mean sojourn time
## G(t, u, t0, mean_sojourn_time): cumulative distribution function of sojourn time and mean sojourn time
##            t: screen time
##            u: time interval between initial age to preclinical onset
##            t0: initial time
##            mean_sojourn_time: mean sojourn time
## h(t, s, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): true sensitivity for both cases
##            t: screen time
##            s: sojourn time, time between preclinical onset and clinical onset
##            u: time interval between initial age to preclinical onset
##            t0: initial time
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## h_g(t, s, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): function of true sensitivity and sojourn time for progressive cases
##            t: screen time
##            u: preclinical onset age
##            s: clinical onset age
##            t0: initial time
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## integral_h_g_progressive(u, t, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): integral of f_g for progressive cases
##            u: preclinical onset age
##            t: screen time
##            s: clinical onset age
##            t0: initial time
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## prospective_sens_num(u, t1, t0, onset_sensitivity, clinical_sensitivity, preonset_rate, mean_sojourn_time):
## numerator of prospective sensitivity
##            u: preclinical onset age
##            t1: upper bound of integral
##            t0: lower bound of integral
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            preonset_rate: exponential rate of cancer onset
##            mean_sojourn_time: mean sojourn time
## prospective_sens_denom(u, t1, t0, preonset_rate, mean_sojourn_time):
## denominator of prospective sensitivity
##            t1: upper bound of integral
##            t0: lower bound of integral
##            preonset_rate: exponential rate of cancer onset
##            mean_sojourn_time: mean sojourn time
##################################################

f <- function(u, preonset_rate=NULL, dist="exponential", interval=NA){
  if(dist=="exponential"){
    return(preonset_rate*exp(-preonset_rate*u))
  }
  if(dist=="uniform"){
    return(1/(interval[2]-interval[1]))
  }
}

g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

G <- function(t, u, t0, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  1-exp(-(1/mean_sojourn_time)*(t-tp))
}

# sensitivity function
h <- function(t, s, u, t0, onset_sensitivity, clinical_sensitivity){
  tp <- t0+u # preclinical onset time
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-tp)
  sensitivity_progressive <- ifelse(sensitivity<=clinical_sensitivity, sensitivity, NA)
  # indolent case
  sensitivity_indolent <- onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/s*(t-tp)))
  
  return(c(sensitivity_progressive, sensitivity_indolent))
}

h_g <- function(s, t, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time, indolent){
  if(indolent){
    return(h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)[2]*g(s, mean_sojourn_time))
  }
  return(h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)[1]*g(s, mean_sojourn_time))
}

f_G <- function(u, t, t0, preonset_rate, mean_sojourn_time, dist, interval){
  f(u, preonset_rate, dist, interval)*(1-G(t, u, t0, mean_sojourn_time))
}

integral_h_g <- function(u, t, t0, preonset_rate, onset_sensitivity, clinical_sensitivity, 
                         mean_sojourn_time, dist, interval, indolent){
  tp <- t0+u # preclinical onset time
  if(indolent){
    f(u, preonset_rate, dist, interval)*
      integrate(Vectorize(h_g),
                lower=0, upper=Inf,
                t=t,
                u=u,
                t0=t0,
                onset_sensitivity=onset_sensitivity,
                clinical_sensitivity=clinical_sensitivity,
                mean_sojourn_time=mean_sojourn_time,
                indolent=indolent)$value
  }else{
    f(u, preonset_rate, dist, interval)*
      integrate(Vectorize(h_g),
                lower=t-tp, upper=Inf,
                t=t,
                u=u,
                t0=t0,
                onset_sensitivity=onset_sensitivity,
                clinical_sensitivity=clinical_sensitivity,
                mean_sojourn_time=mean_sojourn_time,
                indolent=indolent)$value
  }
}

prospective_sens_num <- function(t, t0,
                                 preonset_rate,
                                 onset_sensitivity,
                                 clinical_sensitivity,
                                 mean_sojourn_time,
                                 dist,
                                 interval,
                                 indolent_rate){
  prospective_sens_num_integral_progressive <- (1-indolent_rate) *
    integrate(Vectorize(integral_h_g),
              lower=0, upper=t-t0,
              t=t,
              t0=t0,
              preonset_rate=preonset_rate,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time,
              dist = dist, 
              interval = interval,
              indolent=(indolent_rate>0))$value
  prospective_sens_num_integral_indolent <- indolent_rate *
    integrate(Vectorize(integral_h_g),
              lower=0, upper=t-t0,
              t=t,
              t0=t0,
              preonset_rate=preonset_rate,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time,
              dist = dist, 
              interval = interval,
              indolent=(indolent_rate>0))$value
  
  prospective_sens_num_integral_progressive+prospective_sens_num_integral_indolent
}

prospective_sens_denom <- function(t, t0, preonset_rate, mean_sojourn_time, dist, interval, indolent_rate){
  prospective_sens_denom_integral_progressive <- (1-indolent_rate)*
    integrate(Vectorize(f_G),
              lower=0, upper=t-t0,
              t=t, 
              t0=t0,
              preonset_rate=preonset_rate, 
              mean_sojourn_time=mean_sojourn_time,
              dist = dist, 
              interval = interval)$value
  prospective_sens_denom_integral_indolent <- indolent_rate*
    integrate(Vectorize(f),
              lower=0, upper=t-t0,
              preonset_rate=preonset_rate,
              dist = dist, 
              interval = interval)$value
  
  prospective_sens_denom_integral_progressive+prospective_sens_denom_integral_indolent
}


prospective_sens_analyt <- function(start_age,
                                    test_age,
                                    preonset_rate,
                                    mean_sojourn_time,
                                    onset_sensitivity,
                                    clinical_sensitivity,
                                    dist,
                                    interval,
                                    indolent_rate,
                                    confirmation_test_rate=1,
                                    confirmation_test_sensitivity=1){
  
  numerator <- prospective_sens_num(t=test_age, 
                                    t0=start_age,
                                    preonset_rate=preonset_rate,
                                    onset_sensitivity=onset_sensitivity, 
                                    clinical_sensitivity=clinical_sensitivity, 
                                    mean_sojourn_time=mean_sojourn_time,
                                    dist = dist,
                                    interval = interval,
                                    indolent_rate=indolent_rate)
  
  denomerator <- prospective_sens_denom(t=test_age, 
                                        t0=start_age, 
                                        preonset_rate=preonset_rate, 
                                        mean_sojourn_time=mean_sojourn_time,
                                        dist = dist, 
                                        interval = interval, 
                                        indolent_rate=indolent_rate)    
  
  return(confirmation_test_rate*confirmation_test_sensitivity*numerator/denomerator)
}
