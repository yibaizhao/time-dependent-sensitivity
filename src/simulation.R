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
# datestamp <- "2023-02-27"
# datestamp <- "2023-03-01"
datestamp <- "2023-03-06"

set.seed(1234)

##################################################
# Prospective sensitivity evaluated at specified times
# sojourn_time: sojourn time in years
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
##################################################
test_sensitivity <- function(sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             time=seq(0, 20)){
  alpha=log(1/onset_sensitivity-1)
  beta=1/sojourn_time*(log(1/clinical_sensitivity-1)-alpha)
  tset <- tibble(time=time, sensitivity=1/(1+exp(alpha+beta*time)))
  tset <- tset %>% 
    mutate(sensitivity=case_when(sensitivity>clinical_sensitivity~clinical_sensitivity,
                                 TRUE~sensitivity))
  return(tset)
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
# Simulate natural histories and determine sensitivity
# of a test at the specified age that would be estimated
# in a prospective setting among persons with preclinical
# cancer
# N: number of simulated individuals
# start_age: initial age at start of simulation
# preonset_rate: exponential rate of cancer onset
# mean_sojourn_time: exponential mean sojourn time
# indolent_rate: proportion of non-progressive cancers
# test_age: age at one-time screening test
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
##################################################
prospective_test_sensitivity <- function(N=1000,
                                         start_age=40,
                                         preonset_rate=0.1,
                                         mean_sojourn_time=2,
                                         indolent_rate=0,
                                         test_age=50,
                                         onset_sensitivity=0.2,
                                         clinical_sensitivity=0.8,
                                         indolent_sensitivity=NULL){
  # simulate ages at preclinical onset and clinical diagnosis
  dset <- tibble(onset_age=start_age+rexp(N, rate=preonset_rate),
                 sojourn_time=rexp(N, rate=1/mean_sojourn_time),
                 indolent=rbinom(N, 1, indolent_rate),
                 clinical_age=ifelse(indolent, Inf, onset_age+sojourn_time))
  # confirm retrospective test sensitivity at clinical diagnosis
  cset <- dset %>% filter(is.finite(clinical_age))
  cset <- cset %>% mutate(retro_sens=test_sensitivity(sojourn_time=sojourn_time,
                                                      onset_sensitivity,
                                                      clinical_sensitivity,
                                                      time=sojourn_time)$sensitivity)
  ccheck <- cset %>% with(all(retro_sens == clinical_sensitivity))
  ccheck %>% stopifnot()
  # calculate true prospective test sensitivity among preclinical at test age
  pset <- dset %>% filter(onset_age <= test_age & test_age < clinical_age)
  pset <- pset %>% mutate(prosp_sens=test_sensitivity(sojourn_time=sojourn_time,
                                                      onset_sensitivity,
                                                      clinical_sensitivity,
                                                      time=test_age-onset_age)$sensitivity)
  # update indolent sensitivity
  if(!is.null(indolent_sensitivity)){
    pset <- pset %>% mutate(prosp_sens=case_when(indolent==1~indolent_sensitivity,
                                                 TRUE~prosp_sens))
  }
  pcheck <- pset %>% filter(indolent == 0)
  pcheck <- pcheck %>% with(all(between(prosp_sens,
                                        onset_sensitivity,
                                        clinical_sensitivity)))
  pcheck %>% stopifnot()
  return(list(pset))
}

##################################################
# Visualize bias in true sensitivity obtained in
# retrospective vs prospective screening settings
# dset: tibble of prospective screening settings
# ext: figure filename extension
# saveit: logical indicator of whether to save plot
##################################################
plot_prospective_sensitivity <- function(dset,
                                         clinical_sensitivity,
                                         ext='png',
                                         saveit=FALSE){
  varname <- names(dset)[1]
  varlabel <- switch(varname,
                     test_age='Age at screening test (years)',
                     mean_sojourn_time='Mean sojourn time (years)',
                     indolent_rate='Proportion non-progressive')
  dset <- dset %>% mutate(!!sym(varname):=factor(!!sym(varname)))
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  gg <- ggplot(dset)
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  gg <- gg+geom_boxplot(aes_string(x=varname, y='prosp_sens'),
                        fill='gray',
                        width=1/2,
                        position=position_dodge(width=0.4))
  gg <- gg+scale_x_discrete(name=varlabel)
  gg <- gg+scale_y_continuous(name='Test sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  print(gg)
  if(saveit){
    filname <- str_glue('{varname}_{datestamp}.{ext}')
    ggsave(here("plot", filname),
           gg,
           width=5,
           height=5)
  }
}

##################################################
# Control analysis
##################################################
control <- function(N=10000,
                    preonset_rate=0.1,
                    mean_sojourn_time=c(2, 5, 10),
                    indolent_rate=c(0, 0.2, 0.4),
                    start_age=40,
                    test_age=seq(50, 70, by=10),
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    indolent_sensitivity=0.2,
                    ext='png',
                    saveit=FALSE){
  # visualize how test sensitivity increases over sojourn time
  plot_test_sensitivity(sojourn_time=seq(2, 20, by=2),
                        onset_sensitivity,
                        clinical_sensitivity,
                        ext=ext,
                        saveit=saveit)
  # simulate natural histories
  dset <- expand_grid(test_age, mean_sojourn_time, indolent_rate)
  dset <- dset %>% group_by(test_age, mean_sojourn_time, indolent_rate)
  dset <- dset %>% mutate(results=prospective_test_sensitivity(N=N,
                                                               start_age,
                                                               preonset_rate,
                                                               mean_sojourn_time,
                                                               indolent_rate,
                                                               test_age,
                                                               onset_sensitivity,
                                                               clinical_sensitivity,
                                                               indolent_sensitivity))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  # how does bias depend on testing age (mst 2, 0% non-progressive)?
  dset_test_age <- dset %>% slice_min(mean_sojourn_time)
  dset_test_age <- dset_test_age %>% slice_min(indolent_rate)
  dset_test_age <- dset_test_age %>% select(-mean_sojourn_time, -indolent_rate)
  plot_prospective_sensitivity(dset_test_age,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on mean sojourn time (0% non-progressive, age 50)?
  dset_sojourn_time <- dset %>% slice_min(test_age)
  dset_sojourn_time <- dset_sojourn_time %>% slice_min(indolent_rate)
  dset_sojourn_time <- dset_sojourn_time %>% select(-test_age, -sojourn_time)
  plot_prospective_sensitivity(dset_sojourn_time,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on fraction that are non-progressive (mst 2, age 50)?
  dset_indolent_rate <- dset %>% slice_min(test_age)
  dset_indolent_rate <- dset_indolent_rate %>% slice_min(mean_sojourn_time)
  dset_indolent_rate <- dset_indolent_rate %>% select(-test_age, -mean_sojourn_time)
  plot_prospective_sensitivity(dset_indolent_rate,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
}
control(saveit=TRUE)

