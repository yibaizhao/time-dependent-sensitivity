###########
# Simulate true sensitivity baseline on the function of sojourn time
# Author: Yibai Zhao
###########
library(dplyr)
library(ggplot2)
library(here)

##################
# U: time to preclinical onset, U~f(u)
# Y: sojourn time from preclinical to clinical, Y~g(y)
# S: sensitivity depending on sojourn time, S~h(y)
# test_time: test time
##################

set.seed(234)
N = 1000
rateexample=matrix(c(-.2,.2,0,0,-.16,.16,0,0,0),byrow=T,nrow=3)
rate_mx = rateexample
prop_test_pos <- function(N, rate_mx, test_time, spec=1){
  U <- rexp(N, rate = rate_mx[1,2])
  Y <- rexp(N, rate = rate_mx[2,3])
  U_Y <- U+Y
  dset <- data.frame(ID=1:N,U, Y, U_Y)
  if(test_time){
    # exclude those already clinical at test time and those hasn't been onset
    dset <- dset %>% filter(U_Y>test_time)
    # add true test results
    dset <- dset %>% mutate(true_results = 1*(U<=test_time))
    # estimate sensitivity at sojourn time
    dset <- dset %>% mutate(sens = sens_timeF(test_time-U))
    # get imperfect test results based on sensitivity
    dset <- dset %>% rowwise() %>%
      mutate(sens_results = ifelse(true_results==1, 
                                   rbinom(1, 1, prob = sens),
                                   rbinom(1, 1, prob = 1-spec)))
    # estimate observed sensitivity at test point
    tout <- dset %>% 
      group_by(true_results, sens_results)%>%
      summarize(n=n())%>%
      mutate(prop=n/sum(n))
    
    message("Observed sensitivity at ", test_time, ": ", tout$prop[tout$true_results==1 & tout$sens_results==1])
    return(tout$prop[tout$true_results==1 & tout$sens_results==1])
  }
  # estimate observed sensitivity at clinical onset time
  ## estimate sensitivity at maximum sojourn time
  dset <- dset %>% mutate(sens_clinic_onset = sens_timeF(Y))
  ## get imperfect test results based on sensitivity at maximum sojourn time
  dset <- dset %>% rowwise() %>%
    mutate(sens_clinic_onset_results = rbinom(1, 1, sens_clinic_onset))

  message("Observed sensitivity at clinical onset: ", mean(dset$sens_clinic_onset_results))
  message("Avg of sensitivity at clinical onset: ", mean(dset$sens_clinic_onset))
}
##################
# Function of time dependent sensitivity
# exp(a+bx)/(1+exp(a+bx))
##################
sens_timeF <- function(t){
  1/(1 + exp(-(t)))
}

##################################################
# Control analysis
# Note: method argument controls biomarker estimation
##################################################
control <- function(N,
                    spec=1,
                    saveit){
  rateexample=matrix(c(-.2,.2,0,0,-.16,.16,0,0,0),byrow=T,nrow=3)
  seq_test_time <- seq(1, 10, by = 1)
  out_sens <- sapply(seq_test_time, function(t){prop_test_pos(N=N, rate_mx=rateexample, test_time=t, spec=spec)})
  out <- data.frame(test_time = seq_test_time, sensitivity = out_sens) 
  gg <- out %>% ggplot(aes(x = test_time, y = sensitivity)) +
    geom_point() + geom_smooth(method='loess') +
    xlab("Test time") + ylab("Sensitivity")
  prop_test_pos(N=N, rate_mx=rateexample, test_time=F, spec=spec)
  if(saveit){
    ggsave(here("plot", "sim_sens.png"), gg)
  }
}
control(N=10,
        spec=1,
        saveit=F)


