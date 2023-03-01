##################################################
# Evaluate bias in estimates of a cancer screening
# test obtained in a retrospective setting compared
# to estimates obtained in a prospective setting.
# Basic setup:
# 1. Simulate a population of healthy persons at age 40 years
# 2. Simulate onset (eg, exponential with mean 10 years)
# 3. Specify fraction non-progressive (eg, 0%, 20%, 40%)
# 4. Simulate clinical diagnosis (eg, exponential with mean 2, 5, 10 years)
# 5. Specify a test with sensitivity that increases over sojourn time
#    (e.g., linear from 20% at onset to 80% at clinical diagnosis)
# 6. Simulate one-time screening tests (eg, at ages 50, 60, 70 years)
# 7. Calculate test sensitivity in each prospective screening setting
# Main analyses:
# 1. Confirm that sensitivity estimated in retrospective setting is "correct"
# 2. Confirm that sensitivity estimated in prospective setting is "correct"
# 3. Visualize sensitivity estimated in prospective settings:
#       3a. How does bias depend on testing age (0% non-progressive)?
#       3b. How does bias depend on mean sojourn time (0% non-progressive)?
#       3c. How does bias depend on fraction that are non-progressive?
##################################################
library(tidyverse)
library(purrrlyr)
library(here)
library(foreach)
library(doParallel)

# datestamp <- "2022-12-12"
# datestamp <- "2023-01-26"
# datestamp <- "2023-01-31"
# datestamp <- "2023-02-02"
# datestamp <- "2023-02-06"
# datestamp <- "2023-02-07"
# datestamp <- "2023-02-08"
# datestamp <- "2023-02-09"
# datestamp <- "2023-02-21"
datestamp <- "2023-02-27"

##################################################
# U: time to preclinical onset, U~f(u)
# Y: sojourn time from preclinical to clinical, Y~g(y)
# S: sensitivity depending on sojourn time, S~h(y)
# test_time: test time
##################################################
obs_sensF <- function(N,
                      start_age,
                      mean_preonset_rate,
                      mst,
                      indolent_rate=0,
                      test_age=NULL,
                      alpha,
                      beta_t=NULL,
                      beta_st=NULL){
  # generate individual preclinical onset
  U <- start_age + rexp(N, rate=mean_preonset_rate)
  # generate indicator for indolent cancer
  w <- rbinom(N, 1, indolent_rate)
  # sojourn time
  Y <- rexp(N, rate=1/mst)
  # clinical diagnosed time
  U_Y <- U+Y
  # clinical diagnosed time set to Inf if w=1
  U_Y[w==1] <- Inf
  dset0 <- data.frame(ID=1:N, U, Y, U_Y)
  # get cohorts having preclinical onset but not clinical onset
  dset1 <- dset0 %>% filter(U<=test_age & test_age<U_Y)
  # calculate true sensitivity
  ## preclinical sensitivity
  dset <- dset1 %>% mutate(pre_sens_true=sens_timeF(t=test_age-U, st=Y, alpha=alpha, beta_st=beta_st))
  ## clinical sensitivity
  dset2 <- dset0 %>% mutate(c_sens_true=sens_timeF(t=Y, st=Y, alpha=alpha, beta_st=beta_st))
  # get imperfect test results based on preclinical sensitivity
  dset <- dset %>% by_row(function(row){
                            rbinom(1, 1, prob=row$pre_sens_true)
                      }, .collate="rows", .to="pre_results")
  # get imperfect test results at clinical diagnosed time
  dset2 <- dset2 %>% by_row(function(row){
                              rbinom(1, 1, prob=row$c_sens_true)
                      }, .collate="rows", .to="c_results")
  # calculate mean true sensitivity and observed sensitivity at test age
  dset_out <- data.frame(pre_sens_true=mean(dset$pre_sens_true, na.rm=TRUE),
                         pre_sens_obs=mean(dset$pre_results, na.rm=TRUE),
                         c_sens_true=mean(dset2$c_sens_true, na.rm=TRUE),
                         c_sens_obs=mean(dset2$c_results, na.rm=TRUE))
  dset_out$pre_sens_obs_indolent <- dset_out$pre_sens_obs*indolent_rate
  dset_out$pre_sens_obs_progressive <- dset_out$pre_sens_obs*(1-indolent_rate)
  return(dset_out)
}

##################
# Function of time dependent sensitivity
# exp(a+bx)/(1+exp(a+bx))
## st: sojourn time, time from preclicnial to clinical
##################
sens_timeF <- function(t,
                       st=NULL,
                       alpha,
                       beta_t=NULL,
                       beta_st=NULL){
  if(!is.null(beta_st)){
    return(1/(1+exp(-(alpha+(beta_st/st)*t))))
  }
}

replicate_sensF <- function(B,
                            N,
                            start_age,
                            mean_preonset_rate,
                            mst,
                            test_age,
                            indolent_rate=0,
                            alpha,
                            beta_t=NULL,
                            beta_st=NULL){
  registerDoParallel(8)
  foreach (i=1:B, .combine=rbind) %dopar% {
    obs_sensF(N=N,
              start_age=start_age,
              mean_preonset_rate=mean_preonset_rate,
              mst=mst,
              test_age=test_age,
              indolent_rate=indolent_rate,
              alpha=alpha,
              beta_t=beta_t,
              beta_st=beta_st)
  }
}

outF <- function(B,
                 N,
                 start_age,
                 mean_preonset_rate,
                 mst_seq,
                 test_age_seq,
                 indolent_rate_seq,
                 alpha,
                 beta_t=NULL,
                 beta_st=NULL,
                 saveit=FALSE){
  dset_seq <- expand.grid(test_age_seq, mst_seq, indolent_rate_seq)
  names(dset_seq) <- c("test_age", "mst", "indolent_rate")
  dset_seq <- dset_seq %>% by_row(function(row){
                                    replicate_sensF(B,
                                                    N,
                                                    start_age,
                                                    mean_preonset_rate,
                                                    mst=row$mst,
                                                    test_age=row$test_age,
                                                    indolent_rate=row$indolent_rate,
                                                    alpha,
                                                    beta_t,
                                                    beta_st)
                 }, .collate="row")

  dset_seq <- dset_seq %>% pivot_longer(cols=contains("sens"),
                                        names_to="Type",
                                        values_to="Sensitivity")
  # figure
  labs <- paste0("mst=", mst_seq)
  names(labs) <- mst_seq
  gg_sens <- dset_seq %>%
    ggplot(aes(x=factor(test_age), y=Sensitivity, color=Type)) +
    geom_boxplot(width=1/2, position=position_dodge(width=0.4)) +
    facet_wrap(~mst, ncol=2, labeller=labeller(mst=labs)) +
    xlab("Age at test") +
    scale_color_discrete(name="Sensitivity Type",
                         breaks=unique(dset_seq$Type),
                         labels=c("preclinical true",
                                  "preclinical obs",
                                  "clinical true",
                                  "clinical obs",
                                  "preclinical obs indolent",
                                  "preclinical obs progressive"))
    print(gg_sens)
    if(saveit){
      filname <- str_glue('fig_sens_all_{datestamp}.png')
      ggsave(here("plot", filname),
             gg_sens,
             width=8,
             height=6)
    }
}

plot_sens <- function(t_seq,
                      st_seq=NULL,
                      alpha,
                      beta_t=NULL,
                      beta_st=NULL,
                      saveit=FALSE){
  if(is.null(st_seq) & is.null(beta_st)){
    dset_seq <- data.frame(t=t_seq,
                           sens=sens_timeF(t=t_seq,
                                           alpha=alpha,
                                           beta_t=beta_t))
    gg_sens <- dset_seq %>%
      ggplot(aes(x=t, y=sens)) +
      geom_line() +
      ylim(0, 1) +
      labs(x="Test time (in years)",
           y="Sensitivity")
  } else {
    dset_seq <- expand.grid(t_seq, st_seq)
    names(dset_seq) <- c("t", "st")
    dset_seq <- dset_seq %>% by_row(function(row){
                                      sens_timeF(t=row$t,
                                                 st=row$st,
                                                 alpha=alpha,
                                                 beta_t=beta_t,
                                                 beta_st=beta_st)
           }, .collate="rows", .to="sens")
    gg_sens <- dset_seq %>%
      ggplot(aes(x=t, y=sens, color=st, group=st)) +
      geom_line() +
      ylim(0, 1) +
      labs(x="Test time (in years)",
           y="Sensitivity",
           color='Sojourn time')
  }
  print(gg_sens)
  if(saveit){
    filname <- str_glue('fig_sens_{datestamp}.png')
    ggsave(here("plot", filname), gg_sens)
  }
}

generate_tb <- function(N,
                        mean_st_seq,
                        test_time_seq,
                        beta_t=NULL,
                        beta_st=NULL,
                        saveit=FALSE){
  # estimate observed sensitivity
  out_sens <- sapply(mean_st_seq,
                     function(s) {sapply(test_time_seq,
                                         function(t){obs_sensF(N=N,
                                                               mean_preonset_rate=mean_preonset_rate,
                                                               mean_st=s,
                                                               test_time=t,
                                                               beta_t=beta_t,
                                                               beta_st=beta_st)})})
  out <- data.frame(test_time=test_time_seq, out_sens)
  names(out)[2:(length(mean_st_seq)+1)] <- paste0("mean_st=", mean_st_seq)
  out <- rbind(out, colMeans(out, na.rm=TRUE))
  if(saveit){
    filname <- str_glue('tbl_sens_{datestamp}.csv')
    write.csv(out, here("tables", filname))
  }
  return(out)
}

##################################################
# Control analysis
# Note: method argument controls biomarker estimation
##################################################
control <- function(N=1000,
                    B=10000,
                    mean_preonset_rate=0.1,
                    mean_st=1/0.16,
                    mst_seq=c(2, 5, 10),
                    indolent_rate_seq=c(0, 0.2, 0.4),
                    start_age=40,
                    test_age_seq=seq(50, 70, by=10),
                    alpha=-2.3858,
                    beta_st=3.7721,
                    saveit=FALSE){
  set.seed(234)
  # how does bias depend on testing age (0% non-progressive)?
  # how does bias depend on mean sojourn time (0% non-progressive)?
  # how does bias depend on fraction that are non-progressive?
  plot_sens(t_seq=seq(0, 30, 1),
            st_seq=seq(2, 10, 2),
            alpha=alpha,
            beta_st=beta_st,
            saveit=saveit)
  # compare between preclinical and clinical sensitivity under indolent cases
  outF(B=B,
       N=1000,
       start_age=start_age,
       mean_preonset_rate=mean_preonset_rate,
       mst_seq=mst_seq,
       test_age_seq=test_age_seq,
       indolent_rate_seq=indolent_rate_seq,
       alpha=alpha,
       beta_t=NULL,
       beta_st=beta_st,
       saveit=saveit)
}
control(saveit=TRUE)

