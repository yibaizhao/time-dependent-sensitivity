##################################################
# Compare true sensitivity VS empirical sensitivity in the simulation
##################################################
library(tidyverse)
library(ggplot2)
library(here)
source("simulation_functions_empirical_sensitivity_expanded.R")
# datestamp <- '2023-10-16'
# datestamp <- '2023-10-17'
# datestamp <- '2023-10-24'
# datestamp <- '2023-10-25'
datestamp <- '2023-10-27'

set.seed(1234)

##################################################
# Distribution of true sensitivity
# Progressive cancer: linear regression
# Indolent cancer: a+(b-a)*exp(1-beta*t), where beta is a function of sojourn time
##################################################
test_sensitivity <- function(sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             is_indolent=FALSE,
                             time,
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
      mutate(sensitivity=ifelse(sensitivity<clinical_sensitivity, sensitivity, clinical_sensitivity),
             sensitivity=ifelse(time<0|time>sojourn_time, NA, sensitivity))
  }else{
    tset <- tibble(time=time, 
                   sensitivity=onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/sojourn_time*time)))
  }
  return(tset)
}

##################################################
# Emprirical sensitivity in the simulation
##################################################
# use the breast cancer example, fit to BCSC data
#simulation method
sim_empirical_sensitivity <- function(preonset.rate=0.2, 
                                      mst, 
                                      start.screen.time, 
                                      end.screen.time, 
                                      post.screen.lookout, 
                                      nreps,
                                      multiple.tests){
  rateexample <- matrix(c(-preonset.rate,preonset.rate,0,
                          0,-1/mst,1/mst,
                          0,0,0),
                        byrow=T,nrow=3)
  screen.times <- seq(start.screen.time, end.screen.time, by=post.screen.lookout)
  sens_sim <- empirical.sensitivity.simulation(rate.matrix=rateexample,
                                               start.dist=c(1,0,0),
                                               screen.times=screen.times,
                                               post.screen.lookout=post.screen.lookout,
                                               clinical.cancer.state=3,
                                               pre.clinical.cancer.state=2,
                                               nreps=nreps,
                                               multiple.tests=multiple.tests)
  
  return(sens_sim)
}

plot_empirical_sensitivity <- function(mean.sojourn.time, 
                                       start.screen.time, 
                                       end.screen.time, 
                                       post.screen.lookout, 
                                       nreps,
                                       saveit){
  # get the total number of cores
  numOfCores <- detectCores()
  # register all the cores
  registerDoParallel(numOfCores) 
  oset <- foreach(mst=mean.sojourn.time, .combine=rbind) %dopar%  {
    dset_TRUE <- sim_empirical_sensitivity(preonset.rate=0.2, 
                                      mst, 
                                      start.screen.time, 
                                      end.screen.time, 
                                      post.screen.lookout, 
                                      nreps,
                                      multiple.tests=TRUE)
    dset_TRUE$multiple <- TRUE
    dset_FALSE <- sim_empirical_sensitivity(preonset.rate=0.2, 
                                           mst, 
                                           start.screen.time, 
                                           end.screen.time, 
                                           post.screen.lookout, 
                                           nreps,
                                           multiple.tests=FALSE)
    dset_FALSE$multiple <- FALSE
    dset <- rbind(dset_TRUE, dset_FALSE)
    dset$mst <- mst
    return(dset)
  }
  oset <- oset %>% pivot_longer(c("empirical_sensitivity", "mean_true_sensitivity"),
                                        names_to = "type",
                                        values_to = "sensitivity")
  oset_es <- oset %>% filter(type=="empirical_sensitivity")
  oset_es <- oset_es %>% pivot_longer(c("screen", "interval"),
                                                names_to = "sens_type",
                                                values_to = "count")
  gg <- oset %>%
    ggplot() +
    geom_line(aes(x=obs.times, y=sensitivity, color=type, linetype=multiple)) +
    geom_point(data=oset_es, aes(x=obs.times, y=sensitivity, size=count, color=sens_type), alpha=0.5) +
    facet_grid(.~mst, labeller=label_both) +
    scale_x_continuous(breaks = seq(0, 10, by = 1)) +
    ylim(0, 1)
  print(gg)
  
  if(saveit){
    filename <- str_glue("fig_es_vs_ts_{datestamp}.png")
    ggsave(plot=gg, here("plot", filename))
    
  }
  
}

control <- function(preonset.rate=0.2,
                    mean.sojourn.time=c(1, 2, 4, 6), 
                    start.screen.time=1, 
                    end.screen.time=10, 
                    post.screen.lookout=1, 
                    nreps=10000,
                    saveit){
  plot_empirical_sensitivity(mean.sojourn.time, 
                             start.screen.time, 
                             end.screen.time, 
                             post.screen.lookout, 
                             nreps,
                             saveit)
}

control(saveit=TRUE)

