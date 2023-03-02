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
#    (e.g., logarithmic from 20% at onset to 80% at clinical diagnosis)
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
library(viridis)
library(scales)
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
#datestamp <- "2023-02-27"
datestamp <- "2023-03-01"

##################################################
# Prospective sensitivity evaluated at specified times
# sojourn_time: sojourn time in years
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
##################################################
test_sensitivity <- function(sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             alpha=log(1/onset_sensitivity-1),
                             beta=1/sojourn_time*(log(1/clinical_sensitivity-1)-alpha),
                             time=seq(0, 20)){
  return(tibble(time=time, sensitivity=1/(1+exp(alpha+beta*time))))
}

plot_test_sensitivity <- function(sojourn_time,
                                  onset_sensitivity,
                                  clinical_sensitivity,
                                  ext='png',
                                  saveit=FALSE){
  sset <- tibble(sojourn_time)
  sset <- sset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity))
  sset <- sset %>% unnest(sensitivity)
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  gg <- ggplot(sset)
  gg <- gg+geom_hline(aes(yintercept=onset_sensitivity), linetype='dashed')
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  gg <- gg+geom_line(aes(x=time,
                         y=sensitivity,
                         group=sojourn_time,
                         colour=sojourn_time))
  gg <- gg+scale_x_continuous(name='Years since onset',
                              limits=c(0, max(sojourn_time)),
                              breaks=c(0, sojourn_time),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='Test sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_colour_viridis(name='Sojourn\ntime')
  print(gg)
  if(saveit){
    filname <- str_glue('fig_sens_{datestamp}.{ext}')
    ggsave(here("plot", filname), gg)
  }
}

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
##################################################
control <- function(N=1000,
                    B=10000,
                    mean_preonset_rate=0.1,
                    mean_sojourn_time=c(2, 5, 10),
                    indolent_rate=c(0, 0.2, 0.4),
                    start_age=40,
                    test_age=seq(50, 70, by=10),
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    ext='png',
                    saveit=FALSE){
  set.seed(234)
  # visualize how test sensitivity increases over sojourn time
  plot_test_sensitivity(sojourn_time=seq(2, 20, by=2),
                        onset_sensitivity,
                        clinical_sensitivity,
                        ext=ext,
                        saveit=saveit)
  # simulate natural histories
  dset <- expand_grid(test_age, mean_sojourn_time, indolent_rate)
  stop()
  dset <- dset %>% by_row(function(row){
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
  sset <- sset %>% mutate(sensitivity=map(mean_sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity))
  sset <- sset %>% unnest(sensitivity)
  # how does bias depend on testing age (mst 2, 0% non-progressive)?
  plot_sensitivity(test_age=test_age,
                   mean_sojourn_time=min(mean_sojourn_time),
                   indolent_rate=min(indolent_rate),
                   saveit=saveit)
  # how does bias depend on mean sojourn time (0% non-progressive, age 50)?
  plot_sensitivity(test_age=min(test_age),
                   mean_sojourn_time=mean_sojourn_time,
                   indolent_rate=min(indolent_rate),
                   saveit=saveit)
  # how does bias depend on fraction that are non-progressive (mst 2, age 50)?
  plot_sensitivity(test_age=min(test_age),
                   mean_sojourn_time=min(mean_sojourn_time),
                   indolent_rate=indolent_rate,
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

