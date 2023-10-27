library(msm)
library(dplyr)
library(foreach)
library(doParallel)

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
    numOfCores <- detectCores()
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
    }
    return(obsout %>% filter(obs.times %in% screen.times))
}


###################################################################################################
#Author YZ&JL
# This function ascertains the empirical.sensitivity after a series of screens, using a simulation approach
#empirical sensitvity is defined as the fraction of all cancers that are screeen detected=screen/(screen +interval)
# INPUTS: rate.matrix=transition matrix of cancer model
#         start.dist=starting probability distribution for cancer model
#         sensitivity of the test for detecting pre-clinical cancer
#         screen.times=times of screening
#         post.screen.lookout = duration of time to look out for interval cancer after last screen
#         clinical.cancer.state= state corresponding to clinical cancer
#         pre.clinical.cancer.state=state corresponding to pre-clinical cancer
# OUTPUTS: empirical sensitivity
#
#################################################################################################

empirical.sensitivity.simulation <- function(seed=123, 
                                             rate.matrix,
                                             start.dist,
                                             screen.times,
                                             post.screen.lookout,
                                             clinical.cancer.state,
                                             pre.clinical.cancer.state,
                                             nreps,
                                             multiple.tests){
  
  # out=replicate(n=nreps,get.dx(sensitivity=sensitivity,rate.matrix=rate.matrix,start.dist=start.dist,screen.times=screen.times,
  #                clinical.cancer.state=clinical.cancer.state,
  #               pre.clinical.cancer.state=pre.clinical.cancer.state,post.screen.lookout),simplify=T)
  set.seed(seed)
  # get the total number of cores
  numOfCores <- detectCores()
  # register all the cores
  registerDoParallel(numOfCores) 
  out <- foreach(id=1:nreps, .combine=rbind) %dopar%  {
    dset <- get.dx(rate.matrix=rate.matrix,
                   start.dist=start.dist,
                   screen.times=screen.times,
                   clinical.cancer.state=clinical.cancer.state,
                   pre.clinical.cancer.state=pre.clinical.cancer.state,
                   post.screen.lookout=post.screen.lookout,
                   multiple.tests=multiple.tests)
    dset$id <- id
    dset
  }

  out_sens <- out %>%
    group_by(obs.times) %>%
    summarise(screen=sum(result==2, na.rm=TRUE),
              interval=sum(result==3, na.rm=TRUE),
              empirical_sensitivity=screen/(screen+interval),
              mean_true_sensitivity=mean(true.sens, na.rm = T))
  return(out_sens)
}
