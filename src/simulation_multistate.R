##################################################
# Mathematical formula of empirical sensitivity
# Author: Yibai Zhao
##################################################
library(msm)
library(dplyr)
library(foreach)
library(doParallel)
library(tidyverse)
library(here)
library(readr)
library(ggplot2)

datastamp <- '02-13-2024'

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

get.observed.datapoint<-function(underlying.state,emission.matrix){
  ###################################################################################################
  #Author: 
  # This function gets an observed data point in a HMM based on an underlying state an emission matrix
  #INPUTS: underlying.state = unobserved underlying state in HMM,
  #        emmision.matrix=a matrix with the emission probablities.
  #        the ith row corresponds to the hidden value X(t)=i, and the kth column to O(t)=k|X(t)=i
  #        thus the rows sum to 1, and k columns correspond to the k possible observed states
  #OUTPUTS: the observed data point
  #################################################################################################
  states<-seq(1:dim(emission.matrix)[2])
  probs<-emission.matrix[underlying.state,]
  # browser()
  sample(x=states,size=1,prob=probs)
}


observed.data.hmm<-function(obs.times,underlying.states,emission.matrix){
  ###################################################################################################
  #Author JL
  # This function the observed states in an HMM for multiple observation times
  # INPUTS: obs.times=a vector of observatio11n times; underlying states=corresponding underlying states
  #         emission.matrix=emission matrix for observed data
  #These functins get the observed data point based on an underlying state and a emission matrix
  #
  #
  #################################################################################################
  obs.data<-sapply(underlying.states,FUN="get.observed.datapoint",emission.matrix)
  # browser()
  out<-data.frame(obs.times, obs.data)
  colnames(out)<-c("obs.times","obs.data")
  return(out)
}


get.dx=function(rate.matrix,
                start.dist,
                screen.times,
                post.screen.lookout,
                clinical.cancer.state,
                pre.clinical.cancer.state,
                multiple.tests=TRUE){
  ###################################################################################################
  #Author YZ&JL
  # This function ascertains whether a person is screen or interval detected during a series of screens
  # INPUTS: rate.matrix=transition matrix of cancer model
  #         start.dist=starting probability distribution for cancer model
  #         sensitivity of the test for detecting pre-clinical cancer
  #         screen.times=times of screening
  #         post.screen.lookout = duration of time to look out for interval cancer after last screen
  #         clinical.cancer.state= state corresponding to clinical cancer
  #         pre.clinical.cancer.state=state corresponding to pre-clinical cancer
  # OUTPUTS: result=1 if no cancer discovered, 2 if screen detected, 3=interval detected
  #
  #################################################################################################
  #  browser()
  nstates=dim(rate.matrix)[1]
  thetimes=c(screen.times,max(screen.times)+post.screen.lookout)
  
  # generate true natural history
  the.start.state=sample(1:nstates,size=1,prob=start.dist)
  outctmc=sim.ctmc(start.state=the.start.state,rate.matrix=rate.matrix, 
                   end.time=1000,start.time=0,absorbing.state=clinical.cancer.state)
  sojourn_time <- outctmc$times[clinical.cancer.state]-outctmc$times[pre.clinical.cancer.state]
  
  # generate the true state at observe time
  discreteout=discrete.ctmc(ctmc.times = outctmc$times,ctmc.states=outctmc$states,obs.times=thetimes)
  
  # generate observed state at observe time
  # get the total number of cores
  numOfCores <- detectCores() #- 2
  # register all the cores
  registerDoParallel(numOfCores) 
  obsout <- foreach(t=discreteout$obs.times, .combine=rbind) %dopar%  {
    t.state=discreteout$states[discreteout$obs.times==t]
    sensitivity <- test_sensitivity(sojourn_time=sojourn_time,
                                    onset_sensitivity=0.2,
                                    clinical_sensitivity=0.8,
                                    is_indolent=FALSE,
                                    time=t-outctmc$times[pre.clinical.cancer.state], # test time-preclinical onset time
                                    method='linear')$sensitivity
    # print(paste("t=", t, ", state=", t.state, ", sens=", sensitivity))
    emission.mat=diag(x=1,nrow=nstates,ncol=nstates)
    emission.mat[pre.clinical.cancer.state,pre.clinical.cancer.state]=sensitivity
    emission.mat[pre.clinical.cancer.state,1]=1-sensitivity
    dset_obs <- observed.data.hmm(obs.times=t,underlying.states=t.state,emission.matrix = emission.mat)
    dset_obs$true.sens <- ifelse(sensitivity==0, NA, sensitivity)
    dset_obs
  }
  obsout$state <- discreteout$states
  
  # obtain 1=no cancer detected, 2=screen, 3=interval
  obsout$result <- NA
  for(t in screen.times){
    t_post <- t+post.screen.lookout
    obsout_t <- obsout %>% filter(obs.times %in% c(t, t_post))
    ## state 1 otherwise
    result=1
    ## state 3: if test time is negative, post screen time is state 3; OR test time already state 3
    result[(obsout_t$obs.data[1]==1 & obsout_t$obs.data[2]==clinical.cancer.state) | obsout_t$obs.data[1]==clinical.cancer.state]=3
    ## state 2: if none of above and test time is state 2
    result[obsout_t$obs.data[1]==pre.clinical.cancer.state]=2
    obsout$result[obsout$obs.times==t] <- result
  }
  if(multiple.tests){
    # censored once screen detected or interval cancer
    ## find time first detected or diagnosed
    stop_time <- obsout$obs.times[min(which(obsout$result %in% c(pre.clinical.cancer.state, clinical.cancer.state)))]
    obsout$result[obsout$obs.times>stop_time] <- NA
    obsout$true.sens[obsout$obs.times>stop_time] <- NA
  }else{
    # censored once clinical onset
    ## find time first clinical onset
    stop_time <- obsout$obs.times[min(which(obsout$state == clinical.cancer.state))]
    obsout$result[obsout$obs.times>=stop_time] <- NA
    obsout$true.sens[obsout$obs.times>=stop_time] <- NA
  }
  
  return(obsout %>% filter(obs.times %in% screen.times))
}

# function of true sensitivity at time t
h <- function(t,
              sojourn_time,
              lower,
              upper){
  # t: time from preclincial onset
  beta0 <- lower
  beta1 <- (upper - lower)/sojourn_time
  true_sens <- beta0 + beta1 * t
  true_sens <- ifelse(true_sens<lower | true_sens>upper, NA, true_sens)
  # true_sens <- pmax(lower, pmin(upper, true_sens))
  return(true_sens)
}


sim_stage_prosp_sens <- function(test_time, 
                                 dset,
                                 include_es_ls = FALSE){
  # find state at test time
  dset <- dset %>%
    filter(times < test_time) %>% # Keep rows where value is less than t
    group_by(id) %>%
    slice_max(order_by = times, n = 1) %>% # Select the last row
    ungroup()
  # select in early stage
  dset <- dset %>% filter(states == 2)
  # find true sensitivity at observe time                
  dset <- dset %>%
    mutate(sensitivity = NA,
           # slope of early stage to clinical
           sensitivity = ifelse(es_cl, 
                                h(t = (test_time-preclinical_onset), 
                                  sojourn_time = sojourn_time_slope,
                                  lower = 0,
                                  upper = 0.3),
                                sensitivity), 
           # slope of early stage to late stage
           sensitivity = ifelse(!es_cl, 
                                h(t = (test_time-preclinical_onset), 
                                  sojourn_time = sojourn_time_slope,
                                  lower = 0,
                                  upper = 0.8),
                                sensitivity)
    )
  # fraction of early state to clinical and early state to late stage transition
  dset <- dset %>% mutate(frac_es_cl = mean(es_cl))
  # filename <- str_glue("tbl_m1_{mst_es_cl}_m2_{mst_es_ls}_test_{test_time}_{datastamp}.csv")
  # write_csv(dset, here("tables", filename))
  dset <- read_csv(here("tables", filename))
  
  if(!include_es_ls){
    dset <- dset %>% filter(es_cl)
  }
  return(c(sensitivity = mean(dset$sensitivity),
              fraction = unique(dset$frac_es_cl)))
}

sim_sens <- function(N,
                     test_time,
                     pre_onset_rate,
                     mst_es_cl,
                     mst_es_ls,
                     mst_ls_cl,
                     saveit = TRUE){
  rate.matrix = matrix(c(-pre_onset_rate, pre_onset_rate, 0, 0, 
                         0, -(1/mst_es_cl + 1/mst_es_ls), 1/mst_es_ls, 1/mst_es_cl,
                         0, 0, -1/mst_ls_cl, 1/mst_ls_cl,
                         0, 0, 0, 0),
                       byrow = T, nrow = 4)
  start.dist = c(1, 0, 0, 0)
  early.state = 2
  late.state = 3
  clinical.cancer.state = 4
  
  set.seed(123)
  # get the total number of cores
  numOfCores <- detectCores() #- 2
  # register all the cores
  registerDoParallel(numOfCores) 
  dset_multistate <- foreach(id=1:N, .combine=rbind) %dopar%  {
    # generate true natural history
    outctmc <- sim.ctmc(start.state=1, rate.matrix=rate.matrix, 
                        end.time=1000, start.time=0, absorbing.state=clinical.cancer.state)
    dset <- data.frame(id = id,
                       times = outctmc$times,
                       states = outctmc$states)
    dset <- dset %>% 
      group_by(id) %>%
      mutate(es_cl = all(states != 3)) %>%
      ungroup()
    # find preclinical onset
    dset_pre <- dset %>% 
      group_by(id) %>%
      reframe(preclinical_onset = times[states==2])
    dset <- dset %>% left_join(dset_pre)
    # find sojourn time
    dset <- dset %>% mutate(sojourn_time = (lead(times) - times))
    # find sojourn time in slope, sojourn time of early stage and late stage
    dset_slope <- dset %>% 
      group_by(id) %>%
      mutate(sojourn_time_slope = sojourn_time[states == 2]) %>%
      ungroup()
    dset <- dset %>% left_join(dset_slope)
    
    return(dset)
  }
  
  if(saveit){
    filename <- str_glue("tbl_m1_{mst_es_cl}_m2_{mst_es_ls}_{datastamp}.csv")
    write_csv(dset_multistate, here("tables", filename))
  }
  # prosp_sens_es <- sapply(test_time, function(t) {
  #   sim_stage_prosp_sens(test_time = t, 
  #                        dset = dset_multistate,
  #                        include_es_ls = FALSE)})
  prosp_sens_es_ls <- sapply(test_time, function(t) {
    sim_stage_prosp_sens(test_time = t, 
                         dset = dset_multistate,
                         include_es_ls = TRUE)})
  dset <- data.frame(time = test_time,
                     # prosp_sens_es,
                     sensitivity = prosp_sens_es_ls[1,],
                     fraction = prosp_sens_es_ls[2,])
  
  return(dset)
}

plot_sens <- function(dset, ext = "png", saveit){
  gg <- dset %>%
    filter(mst_es_cl %in% c(2, 6, 10) & mst_es_ls %in% c(2, 6, 10)) %>%
    ggplot(aes(x = time, y = sensitivity, color = factor(mst_es_cl))) +
    geom_point(aes(size = fraction), alpha = 0.5) +
    geom_line() +
    facet_grid(.~mst_es_ls, 
               labeller = labeller(mst_es_ls = label_both)) +
    ylim(0, 0.8) +
    scale_x_continuous(breaks = unique(dset$time)) +
    xlab("Time since study began") +
    ylab("Sensitivity in early stage") +
    scale_color_discrete(name = "MST from early stage to clinical") +
    scale_size_continuous(name = "Fraction of early stage to clinical") +
    theme(legend.position = "bottom",
          legend.box = "vertical")
  print(gg)
  
  if(saveit){
    filename <- str_glue("fig_multistate_{datastamp}.{ext}")
    ggsave(plot = print(gg), here('figure', filename), width = 8, height = 4)
  }
}

control <- function(N = 10000,
                    test_time = 1:10,
                    pre_onset_rate = 0.1,
                    mst_es_cl = seq(2, 10, 2),
                    mst_es_ls = seq(2, 10, 2),
                    mst_ls_cl = 2,
                    saveit){
  set.seed(123)
  cset <- expand_grid(mst_es_cl, mst_es_ls)
  cset$id <- 1:nrow(cset)
  no_cores <- detectCores()# - 2  # leaving one core free
  registerDoParallel(cores=no_cores)
  dset_all <- 
    foreach(i = cset$id, .combine = "rbind") %dopar% {
    dset <- sim_sens(N = N,
             test_time = test_time,
             pre_onset_rate = pre_onset_rate,
             mst_es_cl = cset$mst_es_cl[i],
             mst_es_ls = cset$mst_es_ls[i],
             mst_ls_cl = mst_ls_cl)
    dset <- data.frame(id = i, dset)
    }
  # dset_all <- dset_all %>% right_join(cset)
  # # write_csv(dset_all, here("tables", "tbl_all.csv"))
  dset_all <- read_csv(here("tables", "tbl_all.csv"))
  plot_sens(dset = dset_all, saveit = saveit)
  
  # filename <- str_glue("tbl_m1_6_m2_4_test_5_{datastamp}.csv")
  # dset <- read_csv(here("tables", filename))
  # dset %>% group_by(es_cl) %>% summarise(mean(sensitivity))
}
control(saveit = TRUE)




