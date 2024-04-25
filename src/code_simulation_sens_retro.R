library(tidyverse)
library(foreach)
library(doParallel)

##################################################
# Simulation of retrospective sensitivity
# Author: Yibai Zhao
## t0: initial time
## t_p: pre-clinical onset time
## t_c: clinical onset time
## t2: clinical check time
## f(u): onset time ~ exp(rate)
## g(s): sojourn time ~ exp(1/mst), mst=mean sojourn time
##################################################
# sensitivity function
h <- function(t, s, tp, onset_sensitivity, clinical_sensitivity){
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-tp)
  sensitivity <- ifelse(t<tp | t>(tp+s), NA, sensitivity)
  
  return(sensitivity)
}

########################################################################################
sim.ctmc<-function(start.state,rate.matrix, end.time,start.time=0,absorbing.state=0){
  ########################################################################################
  #Author: JL
  #This function simulates from a homogeneous CTMC characterized by rate.matrix
  #INPUTS: rate.matrix=rate matrix, start.state=starting state for CTMC, end.time=time
  #         to stop data simulations; start.time= time to start data simulations
  #        absorbing.state=possible multiple absorbing states in the model
  #OUTPUTS: a list with two objects: "times"=transition times and "states"=transition states
  #
  #
  ##########################################################################################
  #browser()
  state.space<-seq(1:dim(rate.matrix)[1])
  size<-dim(rate.matrix)[1]
  cur.state<-start.state
  times<-vector()
  states<-vector()
  
  times[1]<-start.time
  states[1]<-cur.state
  cur.time<-start.time
  k<-2
  while((cur.time<end.time)&(!(cur.state%in% absorbing.state))){
    exp.rate<-(-1)*rate.matrix[cur.state,cur.state]
    if(exp.rate==0){
      cur.time<-end.time
    } else{
      cur.time<-cur.time+rexp(n=1,rate=exp.rate)
      if(cur.time<end.time){
        times[k]<-cur.time
        if(size==2){
          cur.state=as.numeric(state.space[-cur.state])
        }else{
          cur.state<-sample(state.space[-cur.state],size=1,prob=rate.matrix[cur.state,-cur.state])
        }
        states[k]<-cur.state
        k<-k+1
        
      }
    }
  }
  return.list<-list(times,states)
  names(return.list)<-c("times","states")
  return(return.list)
}


discrete.ctmc<-function(ctmc.times,ctmc.states,obs.times){
  ###################################################################################################
  #Author: 
  # This function gets the state of a CTMC at different discrete observation times
  #INPUTS: ctmc.times=the transition times for the CTMC
  #        ctmc.states=the states at each of the transition times
  #        obs.times = the discrete observation times
  #
  #OUTPUTS: a dataframe with two columns: obs.times= observation times, states=value of state at obs.times
  #WARNING: if the observation times are outside of max and min transition time, then
  #         the state is assumed to be unchanged from the closest recorded transition time
  #################################################################################################
  out<- data.frame(approx(x=ctmc.times,y=ctmc.states,xout=obs.times,rule=2,f=0,method="constant"))
  colnames(out)<-c("obs.times","states")
  return(out)
}


sim_data_i <- function(id, 
                       start_age,
                       test_age,
                       preonset_rate,
                       mean_sojourn_time,
                       start.dist,
                       onset_sensitivity,
                       clinical_sensitivity,
                       specificity,
                       pre.clinical.cancer.state,
                       clinical.cancer.state){
  
  rate.matrix <- matrix(c(-preonset_rate, preonset_rate, 0,
                          0, -1/mean_sojourn_time, 1/mean_sojourn_time,
                          0, 0, 0),
                        byrow = T,nrow = 3)
  nstates = dim(rate.matrix)[1]
  # generate true natural history
  the.start.state = sample(1:nstates,size = 1,prob = start.dist)
  outctmc = sim.ctmc(start.state = the.start.state,
                   rate.matrix = rate.matrix, 
                   end.time = 1000,
                   start.time = start_age, 
                   absorbing.state = clinical.cancer.state)
  onset_age <- outctmc$times[pre.clinical.cancer.state]
  clinical_age <- outctmc$times[clinical.cancer.state]
  sojourn_time <- clinical_age-onset_age
  # generate the true state at observe time
  discreteout = discrete.ctmc(ctmc.times = outctmc$times,ctmc.states = outctmc$states,obs.times = test_age)
  state <- discreteout$states
  alpha <- h(t = test_age, 
             s = sojourn_time, 
             tp = onset_age, 
             onset_sensitivity, 
             clinical_sensitivity)
  
  result <- ifelse(state == 2, rbinom(1, 1, alpha),
                   ifelse(state == 1, rbinom(1, 1, 1-specificity), NA))
  dset_i <- data.frame(id, start_age, onset_age, clinical_age, sojourn_time, test_age, state=discreteout$states, alpha, result)
  
  return(dset_i)
}

retros_sens_sim <- function(seed, nreps,
                            t_check,
                             start_age,
                             test_age,
                             preonset_rate,
                             mean_sojourn_time,
                             start.dist,
                             onset_sensitivity,
                             clinical_sensitivity,
                             specificity,
                             pre.clinical.cancer.state,
                             clinical.cancer.state){
  set.seed(seed)
  # get the total number of cores
  numOfCores <- detectCores()
  # register all the cores
  registerDoParallel(numOfCores) 
  dset <- foreach(id=1:nreps, .combine = rbind) %dopar%  {
    dset_i <- sim_data_i(id = id, 
                       start_age = start_age,
                       test_age = test_age,
                       preonset_rate = preonset_rate,
                       mean_sojourn_time = mean_sojourn_time,
                       start.dist = start.dist,
                       onset_sensitivity = onset_sensitivity,
                       clinical_sensitivity = clinical_sensitivity,
                       specificity = specificity,
                       pre.clinical.cancer.state = pre.clinical.cancer.state,
                       clinical.cancer.state = clinical.cancer.state)
    dset_i
  }
  
  # find clinical cases before t_check
  dset <- dset %>% mutate(is_clinical_t_check = (clinical_age <= t_check))
  # filter clincial cases only and not clinical yet at screen time
  dset_clinical <- dset %>% filter(is_clinical_t_check & state != 3)
  retros_sens <- mean(dset_clinical$result)
  
  return(retros_sens)
}


retros_sens_sim(seed = 123, 
                    nreps = 1000,
                    t_check = 80,
                    start_age = 40,
                    test_age = 50,
                    preonset_rate = 1/50,
                    mean_sojourn_time = 6,
                    start.dist = c(1, 0, 0),
                    onset_sensitivity = 0,
                    clinical_sensitivity = 0.8,
                    specificity = 1,
                    pre.clinical.cancer.state = 2,
                    clinical.cancer.state = 3)
