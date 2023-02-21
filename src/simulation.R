###########
# Simulate true sensitivity as a function of sojourn time
# Author: Yibai Zhao
###########
library(dplyr)
library(ggplot2)
library(here)
library(purrrlyr)
library(stringr)
library(tidyr)
library(foreach)
library(doParallel)
##################
# U: time to preclinical onset, U~f(u)
# Y: sojourn time from preclinical to clinical, Y~g(y)
# S: sensitivity depending on sojourn time, S~h(y)
# test_time: test time
##################

# datestamp <- "2022-12-12"
# datestamp <- "2023-01-26"
# datestamp <- "2023-01-31"
# datestamp <- "2023-02-02"
# datestamp <- "2023-02-06"
# datestamp <- "2023-02-07"
# datestamp <- "2023-02-08"
# datestamp <- "2023-02-09"
datestamp <- "2023-02-21"

obs_sensF <- function(N, mean_preonset_rate, mean_st, test_time, beta_t=NULL, beta_st=NULL, Model){
  U <- rexp(N, rate = mean_preonset_rate)
  Y <- rexp(N, rate = 1/mean_st)
  U_Y <- U+Y
  dset0 <- data.frame(ID=1:N,U, Y, U_Y)
  if(test_time){
    # exclude those already clinical at test time and those hasn't been onset
    dset <- dset0 %>% filter(U_Y>test_time & test_time>=U)
    # add true test results
    dset <- dset %>% mutate(true_results = 1*(U<=test_time))
    # estimate sensitivity at sojourn time
    if(Model==1){
      dset <- dset %>% mutate(sens = sens_timeF(t = test_time-U, alpha = alpha, beta_t = beta_t, Model = Model))
    }
    if(Model==2){
      dset <- dset %>% mutate(sens = sens_timeF(t = test_time-U, st = Y, alpha = alpha, beta_t = beta_t, beta_st = beta_st, Model = Model))
    }
    if(Model==3){
      dset <- dset %>% mutate(sens = sens_timeF(t = test_time-U, st = Y, alpha = alpha, beta_st = beta_st, Model = Model))
    }
    # get imperfect test results based on sensitivity
    dset <- dset %>% by_row(function(row){
      rbinom(1, 1, prob = row$sens)}, 
      .collate = "rows", .to = "sens_results")
    # estimate observed sensitivity at test point
    message("Observed sensitivity at ", test_time, ": ", mean(dset$sens_results))
    return(mean(dset$sens_results))
  }
  # estimate observed sensitivity at clinical onset time
  ## estimate sensitivity at maximum sojourn time
  if(Model==1){
    dset2 <- dset0 %>% mutate(sens_clinic_onset = sens_timeF(t = Y, alpha = alpha, beta_t = beta_t, Model = Model))
  }
  if(Model==2){
    dset2 <- dset0 %>% mutate(sens_clinic_onset = sens_timeF(t = Y, st = Y, alpha = alpha, beta_t = beta_t, beta_st = beta_st, Model = Model))
  }
  if(Model==3){
    dset2 <- dset0 %>% mutate(sens_clinic_onset = sens_timeF(t = Y, st = Y, alpha = alpha, beta_st = beta_st, Model = Model))
  }
  ## get imperfect test results based on sensitivity at maximum sojourn time
  dset2 <- dset2 %>% rowwise() %>%
    mutate(sens_clinic_onset_results = rbinom(1, 1, sens_clinic_onset))

  message("Observed sensitivity at clinical onset: ", mean(dset2$sens_clinic_onset_results, na.rm = TRUE))
  message("Avg of sensitivity at clinical onset: ", mean(dset2$sens_clinic_onset, na.rm = TRUE))
}

obs_sens_V2F <- function(N, start_age, mean_preonset_rate, mst, 
                         indolent_rate = 0,
                         test_age=NULL, 
                         alpha, beta_t=NULL, beta_st=NULL, Model){
  # generate individual preclinical onset
  U <- start_age + rexp(N, rate = mean_preonset_rate)
  # generate indicator for indolent cancer
  w <- rbinom(N, 1, indolent_rate)
  # sojourn time
  Y <- rexp(N, rate = 1/mst)
  # sojourn time set to Inf if w=1
  Y[w==1] <- Inf
  U_Y <- U+Y
  dset0 <- data.frame(ID=1:N,U, Y, U_Y)
  
  # get cohorts having preclinical onset, but not clinical onset
  dset1 <- dset0 %>% filter(U<=test_age & test_age<U_Y)
  # calculate true sensitivity
  if(Model==3){
    # preclinical sensitivity
    dset <- dset1 %>% mutate(pre_sens_true = sens_timeF(t = test_age-U, st = Y, alpha = alpha, beta_st = beta_st, Model = Model))
    # clinical sensitivity
    dset2 <- dset0 %>% mutate(c_sens_true = sens_timeF(t = Y, st = Y, alpha = alpha, beta_st = beta_st, Model = Model))
  }
  # get imperfect test results based on preclinical sensitivity
  dset <- dset %>% by_row(function(row){
    rbinom(1, 1, prob = row$pre_sens_true)
    }, .collate = "rows", .to = "pre_results")
  # get imperfect test results at clinical diagnosed time
  dset2 <- dset2 %>% by_row(function(row){
    rbinom(1, 1, prob = row$c_sens_true)
  }, .collate = "rows", .to = "c_results")
  # calculate mean true sensitivity and observed sensitivity at test age
  dset_out <- data.frame(pre_sens_true = mean(dset$pre_sens_true, na.rm = T),
                         pre_sens_obs = mean(dset$pre_results, na.rm = T),
                         c_sens_true = mean(dset2$c_sens_true, na.rm = T),
                         c_sens_obs = mean(dset2$c_results, na.rm = T))
  return(dset_out)
  # c(mean(dset$pre_sens_true), mean(dset$pre_results),
  #   mean(dset$c_sens_true), mean(dset$c_results))
}
  
##################
# Function of time dependent sensitivity
# exp(a+bx)/(1+exp(a+bx))
## st: sojourn time, time from preclicnial to clinical
##################
sens_timeF <- function(t, st=NULL, alpha,
                       beta_t=NULL, beta_st=NULL,
                       Model){
  if(Model==1 & !is.null(beta_t)){
    # return( alpha/(alpha + exp(-(beta_t*t))) )
    return( 1/(1 + exp(-(alpha + beta_t*t))) )
  }
  # 1/(1 + exp(-(beta_t*t + beta_st*st)))
  # alpha/(alpha + exp(-(beta_t*t + beta_st*st)))
  if(Model==2 & !is.null(beta_t) & !is.null(beta_st)){
    return(1/(1 + exp(-(alpha + beta_t*t + beta_st*st))))
  }
  if(Model==3  & !is.null(beta_st)){
    return(1/(1 + exp(-(alpha + (beta_st/st)*t))))
  }
}


replicate_sensF <- function(B, N, start_age, 
                            mean_preonset_rate, mst, test_age, 
                            indolent_rate=0,
                            alpha, beta_t=NULL, beta_st=NULL, Model){
  registerDoParallel(8)  
  foreach (i=1:B, .combine=rbind) %dopar% {
    obs_sens_V2F(N = N, start_age = start_age, 
                 mean_preonset_rate = mean_preonset_rate, mst = mst, 
                 test_age = test_age, 
                 indolent_rate = indolent_rate,
                 alpha = alpha, beta_t = beta_t, beta_st = beta_st, 
                 Model = Model)    
  }
}

outF <- function(B, N, start_age, 
                 mean_preonset_rate, mst_seq, test_age_seq, 
                 indolent_rate=0, alpha, beta_t=NULL, beta_st=NULL, Model, saveit){
  dset_seq <- expand.grid(test_age_seq, mst_seq)
  names(dset_seq) <- c("test_age", "mst")
  
  dset_seq2 <- dset_seq %>% 
    by_row(function(row){
      replicate_sensF(B, N, start_age, mean_preonset_rate, 
                      row$mst, row$test_age, indolent_rate, alpha, beta_t, beta_st, Model=3)
      }, .collate = "row")
  
  dset_seq2 <- 
    dset_seq2 %>% pivot_longer(cols = contains("sens"),
                             names_to = "Type",
                             values_to = "Sensitivity") 
  # figure
  labs <- paste0("mst=", mst_seq)
  names(labs) <- mst_seq
  
  gg_sens <- 
    dset_seq2 %>%
    ggplot(aes(x = factor(test_age), y = Sensitivity, color = Type)) +
    geom_boxplot(width = 1/2, position = position_dodge(width = 0.4)) +
    facet_wrap(~mst, ncol = 2, labeller = labeller(mst = labs)) +
    xlab("Age at test") +
    scale_color_discrete(name="Sensitivity Type",
                       breaks= unique(dset_seq2$Type),
                       labels=c("preclinical true", "preclinical obs", "clinical true", "clinical obs"))
  print(gg_sens)
  if(saveit){
    filname <- str_glue('fig_sens_all_m{Model}_{datestamp}.png')
    ggsave(here("plot", filname), gg_sens, width = 8, height = 6)
  }
  
}

plot_sens <- function(t_seq, st_seq=NULL, alpha, beta_t=NULL, beta_st=NULL, Model, saveit){
  if(is.null(st_seq) & is.null(beta_st)){
    dset_seq <- data.frame(t = t_seq,
                       sens = sens_timeF(t=t_seq, alpha=alpha, beta_t=beta_t, Model=Model))
    gg_sens <- dset_seq %>%
      ggplot(aes(x=t, y=sens)) +
      geom_line() +
      ylim(0, 1) +
      labs(x = "Test time (in years)",
           y = "Sensitivity") 
  }else{
    dset_seq <- expand.grid(t_seq, st_seq)
    names(dset_seq) <- c("t", "st")
    dset_seq <- dset_seq %>% 
      by_row(function(row){sens_timeF(t=row$t, st=row$st, alpha=alpha, beta_t=beta_t, beta_st=beta_st, Model=Model)},
             .collate = "rows", .to = "sens")
    gg_sens <- dset_seq %>%
      ggplot(aes(x=t, y=sens, color=st, group=st)) +
      geom_line() +
      ylim(0, 1) +
      labs(x = "Test time (in years)",
           y = "Sensitivity",
           color = 'Sojourn time') 
  }
  print(gg_sens)
  if(saveit){
    filname <- str_glue('fig_sens_m{Model}_{datestamp}.png')
    ggsave(here("plot", filname), gg_sens)
  }
  
}

generate_tb <- function(N, mean_st_seq, test_time_seq, beta_t=NULL, beta_st=NULL, Model, saveit){
  # estimate observed sensitivity
  if(Model==1){
    out_sens <- sapply(mean_st_seq, function(s) {sapply(test_time_seq, function(t){obs_sensF(N=N,
                                                            mean_preonset_rate=mean_preonset_rate, mean_st=s,
                                                            test_time=t, beta_t=beta_t, 
                                                            Model = Model)})})
  }
  if(Model!=1){
    out_sens <- sapply(mean_st_seq, function(s) {sapply(test_time_seq, function(t){obs_sensF(N=N, 
                                                            mean_preonset_rate=mean_preonset_rate, mean_st=s, 
                                                            test_time=t, beta_t=beta_t, beta_st=beta_st, 
                                                            Model = Model)})})
  }
  out <- data.frame(test_time = test_time_seq,
                    out_sens)
  names(out)[2:(length(mean_st_seq)+1)] <- paste0("mean_st=", mean_st_seq)
  out <- rbind(out, colMeans(out, na.rm = T))
  
  if(saveit){
    filname <- str_glue('tbl_sens_m{Model}_{datestamp}.csv')
    write.csv(out, here("tables", filname))
  }
  return(out)
}


##################################################
# Control analysis
# Note: method argument controls biomarker estimation
##################################################
control <- function(N,
                    B,
                    saveit){
  set.seed(234)
  
  mean_preonset_rate <- 0.2
  mean_st <- 1/0.16
  test_time_seq <- seq(0.5, 10, by = 0.5)
  st_seq <- seq(2, 10, 2)
  #############################
  # Scenario I: assuming everyone has same sensitivity distribution (Model 1)
  #############################
  #@t=10, sens=0.8; @t=2, sens=0.2
  # alpha = 0.125
  # betas = 0.3466
  alpha = -1.5322
  beta_t = 0.2919
  # plot sensitivity distribution for different sojourn time
  plot_sens(t_seq = test_time_seq, alpha = alpha, beta_t = beta_t, Model = 1, saveit = saveit)
  # estimate observed sensitivity
  out_sens <- sapply(test_time_seq, function(t){obs_sensF(N=N, 
                                                          mean_preonset_rate=mean_preonset_rate, mean_st=mean_st, 
                                                          test_time=t, 
                                                          beta_t = beta_t, 
                                                          Model = 1)})
  out <- data.frame(test_time = test_time_seq, sensitivity = out_sens) 
  gg <- out %>% ggplot(aes(x = test_time, y = sensitivity)) +
    geom_point() + geom_smooth(method='loess') +
    ylim(0, 1) +
    xlab("Test time") + ylab("Sensitivity")
  if(saveit){
    filename <- str_glue('fig_obs_sens_m1_{datestamp}.png')
    ggsave(here("plot", filename), gg)
  }
  # estimate overall observed sensitivity and sensitivity at clinical diagnosis time
  obs_sensF(N=N, 
            mean_preonset_rate = mean_preonset_rate, mean_st = mean_st, 
            test_time = F, 
            beta_t = beta_t, 
            Model = 1)
  
  
  #############################
  # Scenario II: Assume sensitivity distribution is a function of sojourn time (Model 2)
  #############################
  #@t=10, st=10, sensitivity=0.8
  # alpha = 18.55
  # betas = c(0.3466, -0.5)
  # alpha = 0.4153
  # betas = c(0.2465, -0.1)
  alpha = 0.5523
  beta_t = 0.1234
  beta_st = -0.2
  # plot sensitivity distribution for different sojourn time
  plot_sens(t_seq = test_time_seq, st_seq = st_seq, 
            alpha = alpha,
            beta_t = beta_t, beta_st = beta_st, 
            Model = 2,
            saveit = saveit)
  # estimate observed sensitivity
  out_sens <- sapply(test_time_seq, function(t){obs_sensF(N=N, 
                                                          mean_preonset_rate=mean_preonset_rate, mean_st=mean_st, 
                                                          test_time=t, beta_t = beta_t, beta_st = beta_st,
                                                          Model = 2)})
  out <- data.frame(test_time = test_time_seq, sensitivity = out_sens) 
  gg <- out %>% ggplot(aes(x = test_time, y = sensitivity)) +
    geom_point() + geom_smooth(method='loess') +
    ylim(0, 1) +
    xlab("Test time") + ylab("Sensitivity")
  if(saveit){
    filename <- str_glue('fig_obs_sens_m2_{datestamp}.png')
    ggsave(here("plot", filename), gg)
  }
  # estimate overall observed sensitivity and sensitivity at clinical diagnosis time
  obs_sensF(N=N, 
            mean_preonset_rate = mean_preonset_rate, mean_st = mean_st, 
            test_time = F, 
            beta_t = beta_t, 
            beta_st = beta_st,
            Model = 2)
  
  #############################
  # Scenario III: Assume sensitivity distribution is a function of sojourn time (Model 3)
  #############################
  alpha = -2.3858
  beta_st = 3.7721
  # plot sensitivity distribution for different sojourn time
  plot_sens(t_seq = test_time_seq, st_seq = st_seq, 
            alpha = alpha,
            beta_st = beta_st, 
            Model = 3,
            saveit = saveit)
  # estimate observed sensitivity
  out_sens <- sapply(test_time_seq, function(t){obs_sensF(N=N, 
                                                          mean_preonset_rate=mean_preonset_rate, mean_st=mean_st, 
                                                          test_time=t, 
                                                          beta_t = beta_t, beta_st = beta_st, 
                                                          Model=3)})
  out <- data.frame(test_time = test_time_seq, sensitivity = out_sens) 
  gg <- out %>% ggplot(aes(x = test_time, y = sensitivity)) +
    geom_point() + geom_smooth(method='loess') +
    ylim(0, 1) +
    xlab("Test time") + ylab("Sensitivity")
  if(saveit){
    filename <- str_glue('fig_obs_sens_m3_{datestamp}.png')
    ggsave(here("plot", filename), gg)
  }
  # estimate overall observed sensitivity and sensitivity at clinical diagnosis time
  obs_sensF(N=N, 
            mean_preonset_rate=mean_preonset_rate, mean_st=mean_st, 
            test_time=F, 
            beta_st = beta_st, 
            Model = 3)
  
  #############################
  # Scenario IV: Assume sensitivity distribution is a function of sojourn time (Model 3)
  #############################
  alpha = -2.3858
  beta_st = 3.7721
  # estimate observed sensitivity
  out_sens <- sapply(test_time_seq, function(t){obs_sens_V2F(N = N, start_age = start_age, 
                                                             mean_preonset_rate = mean_preonset_rate, mst = mst, 
                                                             indolent_rate = 0.3,
                                                             test_age = 60, 
                                                             alpha = alpha, beta_st = beta_st, 
                                                             Model = 3)})
  out <- data.frame(test_time = test_time_seq, sensitivity = out_sens) 
  gg <- out %>% ggplot(aes(x = test_time, y = sensitivity)) +
    geom_point() + geom_smooth(method='loess') +
    ylim(0, 1) +
    xlab("Test time") + ylab("Sensitivity")
  #############################
  # Tables of different mean sojourn time
  #############################
  mean_st_seq <- seq(2, 10, by = 2)
  test_time_seq <- seq(0.5, 10, by = 0.5)
  # Model 1
  alpha = -1.5322
  beta_t = 0.2919
  tb1 <- generate_tb(N, mean_st_seq, test_time_seq, beta_t, Model = 1, saveit = saveit)
  # Model 2
  alpha = 0.5523
  beta_t = 0.1234
  beta_st = -0.2
  tb2 <- generate_tb(N, mean_st_seq, test_time_seq, beta_t = beta_t, beta_st = beta_st, Model = 2, saveit = saveit)
  # Model 3
  alpha = -2.3858
  beta_st = 3.7721
  tb3 <- generate_tb(N, mean_st_seq, test_time_seq, beta_st = beta_st, Model = 3, saveit = saveit)
  
  # Comparison between preclinical and clinical sensitivity
  mean_preonset_rate <- 0.2
  test_age_seq <- seq(60, 80, by = 10)
  mst_seq <- seq(2, 10, 2)
  alpha = -2.3858
  beta_st = 3.7721
  outF(B = B, N = 10000, 
       start_age = 50, mean_preonset_rate = mean_preonset_rate, 
       mst_seq, test_age_seq, 
       alpha = alpha, beta_t=NULL, beta_st = beta_st, Model = 3, 
       saveit = F)
  
  # Comparison between preclinical and clinical sensitivity under indolent cases
  mean_preonset_rate <- 0.2
  test_age_seq <- seq(60, 80, by = 10)
  mst_seq <- seq(2, 10, 2)
  alpha = -2.3858
  beta_st = 3.7721
  outF(B = B, N = 10000, 
       start_age = 50, mean_preonset_rate = mean_preonset_rate, 
       mst_seq, test_age_seq, 
       indolent_rate = 0.1,
       alpha = alpha, beta_t=NULL, beta_st = beta_st, Model = 3, 
       saveit = F)
  
}

control(N=1000,
        B=1000,
        saveit=T)


