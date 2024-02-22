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
library(scales)

# datastamp <- '2024-02-13'
# datastamp <- '2024-02-15'
datastamp <- '2024-02-20'

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
                                 include_es = TRUE,
                                 include_ls = TRUE){
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
  # dset <- read_csv(here("tables", filename))
  
  
  if(include_es & !include_ls){
    # include early stage only
    dset <- dset %>% filter(es_cl&include_es)
  }
  if(include_ls & !include_es){
    # include late stage only
    dset <- dset %>% filter(!es_cl&include_ls)
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
                     saveit = FALSE){
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
  prosp_sens_es <- sapply(test_time, function(t) {
    sim_stage_prosp_sens(test_time = t,
                         dset = dset_multistate,
                         include_es = TRUE, 
                         include_ls = FALSE)})
  prosp_sens_ls <- sapply(test_time, function(t) {
    sim_stage_prosp_sens(test_time = t,
                         dset = dset_multistate,
                         include_es = FALSE, 
                         include_ls = TRUE)})
  prosp_sens_es_ls <- sapply(test_time, function(t) {
    sim_stage_prosp_sens(test_time = t, 
                         dset = dset_multistate,
                         include_es = TRUE,
                         include_ls = TRUE)})
  
  dset <- data.frame(time = test_time,
                     sensitivity_es = prosp_sens_es[1,],
                     sensitivity_ls = prosp_sens_ls[1,],
                     sensitivity = prosp_sens_es_ls[1,],
                     fraction = prosp_sens_es_ls[2,])
  
  return(dset)
}

plot_sens <- function(dset, ext = "png", saveit){
  # prepare dset for plot
  pset <- dset %>% filter(mst_es_cl %in% c(2, 6, 10) & mst_es_ls %in% c(2, 6, 10)) 
  pset <- pset %>% mutate(time = time+40)
  # pset <- pset %>% filter(time %in% seq(0, 10, 2)) 
  # categorize fraction
  pset <- pset %>% mutate(frac_gp = case_when(fraction <= 0.2 ~ "Fraction of early stage to clinical onset: <=20%",
                                              fraction > 0.8 ~ "Fraction of early stage to clinical onset: >80%",
                                              TRUE ~ "Fraction of early stage to clinical onset: >20% and <=80%"))
  
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  
  # gg <- pset %>%
  #   ggplot(aes(x = time, y = sensitivity, color = factor(mst_es_cl), linetype = factor(mst_es_ls))) +
  #   geom_line(size = 1) +
  #   facet_grid(.~frac_gp) +
  #   scale_x_continuous(breaks = unique(pset$time)) +
  #   xlab("Age at screen test") +
  #   ylab("Sensitivity in early stage") +
  #   scale_color_discrete(name = "MST from early stage to clinical") +
  #   scale_linetype_manual(values = c("2" = "dotted", "6" = "dashed", "10" = "solid")) +
  #   scale_linetype_discrete(name = "MST from early stage to late stage") +
  #   theme(legend.position = "bottom",
  #         legend.box = "vertical") +
  #   scale_y_continuous(limits=c(0, 0.8),
  #                      breaks=seq(0, 1, by=0.2),
  #                      labels=label_percent(accuracy=1),
  #                      expand=c(0, 0))
  # print(gg)
  
  pset2 <- pset %>% filter(time == 55)
  # write_csv(pset2, here("tables", "tbl_test_55.csv"))
  gg2 <- pset2 %>%
    ggplot(aes(x = mst_es_cl, y = sensitivity, color = factor(mst_es_ls))) +
    geom_point(aes(size = fraction), alpha = 0.8) +
    geom_line() +
    # scale_x_continuous(breaks = unique(pset$time)) +
    xlab("MST from early stage to clinical") +
    ylab("Sensitivity in early stage") +
    scale_color_discrete(name = "MST from early stage to late stage") +
    scale_size_continuous(name = "Fraction of early stage to clinical") +
    ggtitle("Test at 55-year-old") +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    scale_y_continuous(limits=c(0, 0.8),
                       breaks=seq(0, 1, by=0.2),
                       labels=label_percent(accuracy=1),
                       expand=c(0, 0))
  print(gg2)
  
  gg3 <- pset2 %>%
    ggplot(aes(x = mst_es_cl, y = sensitivity, color = factor(mst_es_ls))) +
    geom_ribbon(aes(ymin = sensitivity_es, ymax = sensitivity_ls, fill = factor(mst_es_ls)), alpha = 0.3, show.legend = FALSE) +
    geom_point(aes(size = fraction), alpha = 0.8) +
    geom_line() +
    # scale_x_continuous(breaks = unique(pset$time)) +
    xlab("MST from early stage to clinical") +
    ylab("Sensitivity in early stage") +
    scale_color_discrete(name = "MST from early stage to late stage") +
    scale_size_continuous(name = "Fraction of early stage to clinical") +
    ggtitle("Test at 55-year-old") +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    scale_y_continuous(limits=c(0, 0.8),
                       breaks=seq(0, 1, by=0.2),
                       labels=label_percent(accuracy=1),
                       expand=c(0, 0))
  print(gg3)
  
  
  if(saveit){
    # filename1 <- str_glue("fig_multistate_V1_{datastamp}.{ext}")
    # ggsave(plot = print(gg), here('figure', filename1), width = 8, height = 4)
    filename2 <- str_glue("fig_multistate_{datastamp}.{ext}")
    ggsave(plot = print(gg2), here('figure', filename2), width = 6, height = 6)
    filename3 <- str_glue("fig_multistate_boundary_{datastamp}.{ext}")
    ggsave(plot = print(gg3), here('figure', filename3), width = 6, height = 6)
  }
}

control <- function(N = 10000,
                    test_time = 10:20,
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
  dset_all <- dset_all %>% right_join(cset)
  write_csv(dset_all, here("tables", "tbl_all.csv"))
  # dset_all <- read_csv(here("tables", "tbl_all_2024-02-22.csv"))
  plot_sens(dset = dset_all, saveit = saveit)
  
  # filename <- str_glue("tbl_m1_6_m2_4_test_5_{datastamp}.csv")
  # dset <- read_csv(here("tables", filename))
  # dset %>% group_by(es_cl) %>% summarise(mean(sensitivity))
}
control(saveit = TRUE)




