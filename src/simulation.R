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
library(stringr)
library(here)
library(foreach)
library(doParallel)
library(tidyr)
library(viridis)

# datestamp <- '2022-12-12'
# datestamp <- '2023-01-26'
# datestamp <- '2023-01-31'
# datestamp <- '2023-02-02'
# datestamp <- '2023-02-06'
# datestamp <- '2023-02-07'
# datestamp <- '2023-02-08'
# datestamp <- '2023-02-09'
# datestamp <- '2023-02-21'
# datestamp <- '2023-02-27'
# datestamp <- '2023-03-01'
# datestamp <- '2023-03-06'
# datestamp <- '2023-04-10'
# datestamp <- '2023-04-18'
# datestamp <- '2023-05-24'
# datestamp <- '2023-06-02'
#datestamp <- '2023-06-05'
#datestamp <- '2023-06-07'
# datestamp <- '2023-06-08'
#datestamp <- '2023-06-09'
# datestamp <- '2023-06-12'
datestamp <- '2023-06-21'

set.seed(1234)

##################################################
# Prospective sensitivity evaluated at specified times
# sojourn_time: sojourn time in years
# onset_sensitivity: test sensitivity at onset
# clinical_sensitivity: test sensitivity at clinical diagnosis
# TODO: Determine annotation positions programmatically
##################################################
test_sensitivity <- function(sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             is_indolent=FALSE,
                             param_indolent=5,
                             time=seq(0, 20),
                             method='linear'){
  # alpha=log(1/onset_sensitivity-1)
  # beta=1/sojourn_time*(log(1/clinical_sensitivity-1)-alpha)
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
      mutate(sensitivity=ifelse(sensitivity>clinical_sensitivity, NA, sensitivity))
  }else{
    tset <- onset_sensitivity+(clinical_sensitivity-onset_sensitivity)/(1+exp(-(1/sojourn_time)*(time-param_indolent)))
  }
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
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               legend.position='bottom')
  gg <- ggplot(sset)
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  gg <- gg+geom_hline(aes(yintercept=onset_sensitivity), color='red')
  gg <- gg+annotate('segment',
                    x=5,
                    xend=5,
                    y=0.1,
                    yend=0.18,
                    colour='red',
                    arrow=arrow(length=unit(0.2, 'cm'), type='closed'))
  gg <- gg+annotate('text',
                    x=5,
                    y=0.1,
                    color='red',
                    vjust=1.5,
                    label='Indolent cancer')
  gg <- gg+geom_segment(data=sset %>% filter(sojourn_time == time),
                        aes(x=5,
                            xend=time,
                            y=0.95,
                            yend=0.85),
                        colour='black',
                        arrow=arrow(length=unit(0.2, 'cm'), type='closed'))
  gg <- gg+annotate('text',
                    x=5,
                    y=0.95,
                    color='black',
                    vjust=-0.5,
                    label='Progressive cancer')
  gg <- gg+geom_line(aes(x=time,
                         y=sensitivity,
                         group=sojourn_time,
                         colour=sojourn_time))
  gg <- gg+scale_x_continuous(name='Years since onset',
                              limits=c(0, max(sojourn_time)),
                              breaks=c(0, sojourn_time),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='True sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_colour_viridis(name='Sojourn time')
  print(gg)
  if(saveit){
    filename <- str_glue('true_sensitivity_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=6,
           width=6)
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
prospective_test_sensitivity <- function(N,
                                         start_age,
                                         preonset_rate,
                                         mean_sojourn_time,
                                         indolent_rate,
                                         test_age,
                                         onset_sensitivity,
                                         clinical_sensitivity,
                                         indolent_sensitivity=NULL,
                                         confirmation_test_rate=NULL,
                                         confirmation_test_sensitivity=NULL){
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
  pset <- pset %>% mutate(prosp_sens_all=test_sensitivity(sojourn_time=sojourn_time,
                                                      onset_sensitivity,
                                                      clinical_sensitivity,
                                                      time=test_age-onset_age)$sensitivity)
  # update indolent sensitivity
  if(!is.null(indolent_sensitivity)){
    pset <- pset %>% mutate(prosp_sens_all=case_when(indolent==1~indolent_sensitivity,
                                                     TRUE~prosp_sens_all))
  }else{
    pset <- pset %>% mutate(prosp_sens_all=case_when(indolent==1~test_sensitivity(sojourn_time=sojourn_time,
                                                                                  onset_sensitivity=onset_sensitivity,
                                                                                  clinical_sensitivity=clinical_sensitivity,
                                                                                  is_indolent=TRUE,
                                                                                  time=(test_age-onset_age))$sensitivity,
                                                     TRUE~prosp_sens_all))
  }
  pcheck <- pset %>% filter(indolent == 0)
  # update confirmation test
  if(!is.null(confirmation_test_rate)&!is.null(confirmation_test_sensitivity)){
    if(confirmation_test_rate>0){
      pset <- pset %>% mutate(prosp_sens_all=prosp_sens_all*confirmation_test_rate*confirmation_test_sensitivity)
    }
  }
  # estimate prospective sensitivity for clinical significant cancer
  pset <- pset %>% mutate(prosp_sens_sig=ifelse(indolent==0, (1-indolent_rate)*prosp_sens_all, NA))

  pcheck <- pcheck %>% with(all(between(prosp_sens_all,
                                        onset_sensitivity,
                                        clinical_sensitivity)))
  pcheck %>% stopifnot()
  return(list(pset))
}

##################################################
# Mathematical formula of prospective sensitivity
## f(t, t0, preonset_rate): function of preclinical onset time and onset rate
##            t: preclinical onset age
##            t0: initial age at start of simulation
##            preonset_rate: preclinical onset rate
## g(t, st, mean_sojourn_time): function of sojourn time
##            t: preclinical onset age
##            st: clinical onset age
##            mean_sojourn_time: mean sojourn time
## G(t, mean_sojourn_time): cumulative distribution function of sojourn time and mean sojourn time
##            t: preclinical onset age
##            t1: screen age
##            mean_sojourn_time: mean sojourn time
## h(t, st, t1, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): true sensitivity
##            t: preclinical onset age
##            st: clinical onset age
##            t1: screen age
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## f_sens_st(s, u, t1, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): function of true sensitivity and sojourn time
##            s: clinical onset age
##            u: preclinical onset age
##            t1: screen age
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## prospective_sens_num(u, t1, t0, onset_sensitivity, clinical_sensitivity, preonset_rate, mean_sojourn_time):
## numerator of prospective sensitivity
##            u: preclinical onset age
##            t1: upper bound of integral
##            t0: lower bound of integral
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            preonset_rate: exponential rate of cancer onset
##            mean_sojourn_time: mean sojourn time
## prospective_sens_denom(u, t1, t0, preonset_rate, mean_sojourn_time):
## denominator of prospective sensitivity
##            t1: upper bound of integral
##            t0: lower bound of integral
##            preonset_rate: exponential rate of cancer onset
##            mean_sojourn_time: mean sojourn time
##################################################

f <- function(t, t0, preonset_rate){
  preonset_rate*exp(-preonset_rate*(t-t0))
}

g <- function(t, st, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*(st-t))
}

G <- function(t, t1, mean_sojourn_time){
  1-exp(-(1/mean_sojourn_time)*(t1-t))
}

h <- function(t, st, t1, onset_sensitivity, clinical_sensitivity){
  alpha <- -log(clinical_sensitivity/onset_sensitivity)
  beta <- -log(onset_sensitivity/clinical_sensitivity)-alpha
  sensitivity <- 1/(1+exp(-(alpha+beta/(st-t)*(t1-t))))
  sensitivity <- ifelse(sensitivity>clinical_sensitivity, clinical_sensitivity, sensitivity)
  return(sensitivity)
}

f_sens_mst <- function(s, u, t1, onset_sensitivity, clinical_sensitivity, mean_sojourn_time){
  h(u, s, t1, onset_sensitivity, clinical_sensitivity)*g(u, s, mean_sojourn_time)
}

prospective_sens_num <- function(u,
                                 t1,
                                 t0,
                                 onset_sensitivity,
                                 clinical_sensitivity,
                                 preonset_rate,
                                 mean_sojourn_time){
  f(u, t0, preonset_rate)*
    integrate(f_sens_mst,
              lower=t1, upper=Inf,
              u=u,
              t1=t1,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time)$value
}

prospective_sens_denom <- function(u, t1, t0, preonset_rate, mean_sojourn_time){
  f(u, t0, preonset_rate)*(1-G(u, t1, mean_sojourn_time))
}


prospective_sens <- function(start_age,
                             screen_age,
                             onset_sensitivity,
                             clinical_sensitivity,
                             preonset_rate,
                             mean_sojourn_time){
  prospective_sens_num_integral <-
    integrate(Vectorize(prospective_sens_num),
              lower=start_age, upper=screen_age,
              t1=screen_age,
              t0=start_age,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              preonset_rate=preonset_rate,
              mean_sojourn_time=mean_sojourn_time)$value
  prospective_sens_denom_integral <-
    integrate(prospective_sens_denom,
              lower=start_age,
              upper=screen_age,
              t1=screen_age,
              t0=start_age,
              preonset_rate=preonset_rate, mean_sojourn_time=mean_sojourn_time)$value
  return(ifelse(prospective_sens_denom_integral==0,
                onset_sensitivity,
                prospective_sens_num_integral/prospective_sens_denom_integral))
}

##################################################
# Compare prospective sensitivity estimated from simulation,
# analytic formula, with true sensitivity
## N: sample size
## start_age: initial age at start of simulation
## test_age: age at one-time screening test
## preonset_rate: exponential rate of cancer onset
## mean_sojourn_time: exponential mean sojourn time
## onset_sensitivity: test sensitivity at onset
## clinical_sensitivity: test sensitivity at clinical diagnosis
##################################################
sens_compare <- function(N,
                         start_age=40,
                         screen_age,
                         preonset_rate=0.1,
                         mean_sojourn_time,
                         onset_sensitivity=0.2,
                         clinical_sensitivity=0.8){
  analytic_sens <- prospective_sens(start_age,
                                    screen_age,
                                    onset_sensitivity,
                                    clinical_sensitivity,
                                    preonset_rate,
                                    mean_sojourn_time)
  sim_sens_tb <- prospective_test_sensitivity(N,
                                              start_age,
                                              preonset_rate,
                                              mean_sojourn_time,
                                              indolent_rate=0,
                                              test_age=screen_age,
                                              onset_sensitivity,
                                              clinical_sensitivity,
                                              indolent_sensitivity=NULL)[[1]]
  sim_sens <- mean(sim_sens_tb$prosp_sens_all)
  return(list(tibble(analytic_sens, sim_sens)))
}

##################################################
# Visualize comparison of prospective sensitivity
# estimated from simulation and analytic formula
## N: sample size
## test_age: age at one-time screening test
## mean_sojourn_time: exponential mean sojourn time
## ext: figure filename extension
## saveit: logical indicator of whether to save plot
##################################################
plot_sens_compare <- function(N,
                              test_age,
                              mean_sojourn_time,
                              ext='png',
                              saveit=FALSE){
  dset <- expand_grid(test_age, mean_sojourn_time)
  dset <- dset %>% group_by(test_age, mean_sojourn_time)
  dset <- dset %>% mutate(results=sens_compare(N=N,
                                               screen_age=test_age,
                                               mean_sojourn_time=mean_sojourn_time))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  dset <- dset %>% pivot_longer(cols=ends_with('sens'),
                                names_to='type', values_to='sensitivity')
  # dset$test_age <- factor(dset$test_age)
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  gg <- ggplot(dset, aes(x=test_age, y=sensitivity, color=type))
  # gg <- gg+geom_boxplot()
  gg <- gg+geom_line()
  gg <- gg+facet_grid(.~mean_sojourn_time, labeller=label_both)
  gg <- gg+xlab('Test age')
  gg <- gg+ylab('Sensitivity')
  gg <- gg+ylim(0,1)
  gg <- gg+guides(color=guide_legend(title='Types of prospective sensitivity'))
  gg <- gg+scale_color_discrete(labels=c('Analytic', 'Simulation'))
  gg <- gg+theme(legend.position='bottom')
  print(gg)

  if(saveit){
    filename <- str_glue('fig_sens_compare_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           width=6,
           height=5)
  }
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
                                         significant_cancer_only=FALSE,
                                         ext='png',
                                         saveit=FALSE){
  varname <- names(dset)[1]
  varlabel <- switch(varname,
                     test_age='Age at screening test (years)',
                     mean_sojourn_time='Mean sojourn time (years)',
                     indolent_rate='Proportion non-progressive',
                     confirmation_test_rate='Percent of subjects who took confirmation test',
                     confirmation_test_sensitivity='Confirmation test sensitivity')
  dset <- dset %>% mutate(!!sym(varname):=factor(!!sym(varname)))
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  gg <- ggplot(dset)
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  if(significant_cancer_only){
    gg <- gg+geom_boxplot(aes_string(x=varname, y='prosp_sens_sig'),
                          fill='gray',
                          width=1/2,
                          position=position_dodge(width=0.4))
    gg <- gg+ylab('Prospective sensitivity of detecting progressive cancer')
  } else {
    gg <- gg+geom_boxplot(aes_string(x=varname, y='prosp_sens_all'),
                          fill='gray',
                          width=1/2,
                          position=position_dodge(width=0.4))
    gg <- gg+ylab('Prospective sensitivity')
  }
  gg <- gg+scale_x_discrete(name=varlabel)
  gg <- gg+scale_y_continuous(limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  print(gg)
  if(saveit){
    if(significant_cancer_only){
      filename <- str_glue('{varname}_sig_{datestamp}.{ext}')
    } else {
      filename <- str_glue('{varname}_{datestamp}.{ext}')
    }
    ggsave(here('plot', filename),
           plot=gg,
           width=5,
           height=5)
  }
}

plot_prospective_sensitivity_compare <- function(N,
                                                 dset,
                                                 start_age,
                                                 preonset_rate,
                                                 onset_sensitivity,
                                                 clinical_sensitivity,
                                                 indolent_sensitivity,
                                                 ext='png',
                                                 significant_cancer_only=FALSE,
                                                 saveit=FALSE){
  # analytic_sens <- prospective_sens(start_age, screen_age, onset_sensitivity, clinical_sensitivity, preonset_rate, mean_sojourn_time)
  dset <- dset %>% group_by(across(all_of(names(dset))))
  dset <- dset %>% mutate(prosp_sens_all=prospective_test_sensitivity(N=N,
                                                                      start_age=start_age,
                                                                      preonset_rate=preonset_rate,
                                                                      mean_sojourn_time=mean_sojourn_time,
                                                                      indolent_rate=indolent_rate,
                                                                      test_age=test_age,
                                                                      onset_sensitivity=onset_sensitivity,
                                                                      clinical_sensitivity=clinical_sensitivity,
                                                                      indolent_sensitivity=indolent_sensitivity,
                                                                      confirmation_test_rate=confirmation_test_rate,
                                                                      confirmation_test_sensitivity=confirmation_test_sensitivity))
  if(significant_cancer_only){
    dset <- dset %>% mutate(mean_prosp_sens=map(prosp_sens_all, ~mean(.x$prosp_sens_sig, na.rm=TRUE)))
    ylabel <- 'Prospective sensitivity for progressive cancer'
  } else {
    dset <- dset %>% mutate(mean_prosp_sens=map(prosp_sens_all, ~mean(.x$prosp_sens_all)))
    ylabel <- 'Prospective sensitivity for any cancer'
  }
  dset <- dset %>% select(-prosp_sens_all)
  dset <- dset %>% ungroup()
  dset <- dset %>% unnest(mean_prosp_sens)
  dset <- dset %>% mutate(degradation=clinical_sensitivity-mean_prosp_sens,
                          test_age=factor(test_age),
                          mean_sojourn_time=factor(mean_sojourn_time),
                          indolent_rate=paste('Indolent\nfraction', sprintf('%2.0f%%', 100*indolent_rate)),
                          indolent_rate=factor(indolent_rate))
                          # confirmation_test_sensitivity=paste('Confirmation test\nsensitivity', sprintf('%2.0f%%', 100*confirmation_test_sensitivity)))
                          #confirmation_test_sensitivity=factor(confirmation_test_sensitivity,
                          #                                     levels=c('Confirmation test\nsensitivity 100%',
                          #                                              'Confirmation test\nsensitivity 80%',
                          #                                              'Confirmation test\nsensitivity 60%')))
  dset <- dset %>% mutate( 
    confirmation_test_rate=factor(paste('Confirmation test\nfrequency/sensitivity', 
                                        sprintf('%2.0f%%', 100*confirmation_test_rate)),
                                  levels=paste('Confirmation test\nfrequency/sensitivity', 
                                               sprintf('%2.0f%%', 100*sort(unique(dset$confirmation_test_rate), decreasing = TRUE)))
                                  ))
  
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               panel.grid.major.y=element_line(),
               panel.spacing=unit(0.04, 'npc'),
               strip.text.y=element_text(color='black', angle=0),
               strip.background=element_rect(color=NA, fill=NA))
  gg <- dset %>% ggplot()
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  gg <- gg+geom_bar(aes(x=test_age,
                        y=mean_prosp_sens,
                        fill=mean_sojourn_time),
                    stat='identity',
                    position=position_dodge(width=1))
  gg <- gg+geom_linerange(aes(xmin=test_age,
                              ymin=mean_prosp_sens,
                              x=test_age,
                              y=mean_prosp_sens,
                              xmax=test_age,
                              ymax=clinical_sensitivity,
                              color=mean_sojourn_time),
                          alpha=0.3,
                          position=position_dodge(width=1))
  gg <- gg+geom_label(aes(x=test_age,
                          y=mean_prosp_sens+degradation/2,
                          label=sprintf('-%2.0f%%', 100*degradation),
                          color=mean_sojourn_time),
                      size=2,
                      label.size=0,
                      label.padding=unit(0.1, 'lines'),
                      position=position_dodge2(width=1),
                      show.legend=FALSE)
  gg <- gg+ylab(ylabel)
  gg <- gg+geom_hline(aes(yintercept=0))
  if(significant_cancer_only){
    gg <- gg+facet_grid(.~indolent_rate)
  } else {
    gg <- gg+facet_grid(indolent_rate~confirmation_test_rate)
  }
  gg <- gg+scale_x_discrete(name='Age at screening test (years)',
                            expand=c(0, 0))
  gg <- gg+scale_y_continuous(limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_color_viridis(name='Mean sojourn\ntime (years)', discrete=TRUE)
  gg <- gg+scale_fill_viridis(name='Mean sojourn\ntime (years)', discrete=TRUE)
  print(gg)
  if(saveit){
    if(significant_cancer_only){
      filename <- str_glue('plot_compare_sig_setting_{datestamp}.{ext}')
    } else {
      filename <- str_glue('plot_compare_setting_{datestamp}.{ext}')
    }
    ggsave(here('plot', filename),
           plot=gg,
           width=10,
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
                    confirmation_test_rate=c(0.4, 0.6, 1),
                    confirmation_test_sensitivity=c(0.6, 0.8, 1),
                    start_age=40,
                    test_age=seq(50, 70, by=10),
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    indolent_sensitivity=NULL,
                    ext='png',
                    saveit=FALSE){
  # visualize how test sensitivity increases over sojourn time
  plot_test_sensitivity(sojourn_time=seq(2, 10, by=2),
                        onset_sensitivity,
                        clinical_sensitivity,
                        ext=ext,
                        saveit=saveit)
  # simulate natural histories
  dset <- expand_grid(test_age, mean_sojourn_time, indolent_rate)
  dset <- dset %>% group_by(test_age, mean_sojourn_time, indolent_rate)
  dset <- dset %>% mutate(results=prospective_test_sensitivity(N=N,
                                                               start_age=start_age,
                                                               preonset_rate=preonset_rate,
                                                               mean_sojourn_time=mean_sojourn_time,
                                                               indolent_rate=indolent_rate,
                                                               test_age=test_age,
                                                               onset_sensitivity=onset_sensitivity,
                                                               clinical_sensitivity=clinical_sensitivity,
                                                               indolent_sensitivity=indolent_sensitivity))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  # how does bias depend on testing age (mst 5, 0% non-progressive)?
  dset_test_age <- dset %>% filter(mean_sojourn_time==5)
  dset_test_age <- dset_test_age %>% filter(indolent_rate==0)
  dset_test_age <- dset_test_age %>% select(-mean_sojourn_time, -indolent_rate)
  plot_prospective_sensitivity(dset_test_age,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on mean sojourn time (0% non-progressive, age 50)?
  dset_sojourn_time <- dset %>% filter(test_age==50)
  dset_sojourn_time <- dset_sojourn_time %>% filter(indolent_rate==0)
  dset_sojourn_time <- dset_sojourn_time %>% select(-test_age, -sojourn_time)
  plot_prospective_sensitivity(dset_sojourn_time,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on fraction that are non-progressive (mst 2, age 50)?
  dset_indolent_rate <- dset %>% filter(test_age==50)
  dset_indolent_rate <- dset_indolent_rate %>% filter(mean_sojourn_time==2)
  dset_indolent_rate <- dset_indolent_rate %>% select(-test_age, -mean_sojourn_time)
  plot_prospective_sensitivity(dset_indolent_rate,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on fraction that are non-progressive in terms of progressive cancer (mst 2, age 50)?
  plot_prospective_sensitivity(dset_indolent_rate,
                               clinical_sensitivity,
                               significant_cancer_only=TRUE,
                               ext=ext,
                               saveit=saveit)
  # how does bias depend on confirmation test (mst 2, age 50, indolent_rate=0.2, confirmation test sensitivity = 1)?
  dset <- expand_grid(test_age=50,
                      mean_sojourn_time=2,
                      indolent_rate=0,
                      confirmation_test_rate,
                      confirmation_test_sensitivity)
  dset <- dset %>% group_by(confirmation_test_rate, confirmation_test_sensitivity)
  dset <- dset %>% mutate(results=prospective_test_sensitivity(N=N,
                                                               start_age,
                                                               preonset_rate,
                                                               mean_sojourn_time,
                                                               indolent_rate,
                                                               test_age,
                                                               onset_sensitivity,
                                                               clinical_sensitivity,
                                                               indolent_sensitivity,
                                                               confirmation_test_rate,
                                                               confirmation_test_sensitivity))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  ## how does bias depend on confirmation test rate (mst 2, age 50, indolent_rate=0.2, confirmation test sensitivity = 1)?
  dset_c_rate <- dset %>% filter(confirmation_test_sensitivity==1)
  dset_c_rate <- dset_c_rate %>% select(-c(test_age,
                                           mean_sojourn_time,
                                           indolent_rate,
                                           confirmation_test_sensitivity))
  plot_prospective_sensitivity(dset_c_rate,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  ## how does bias depend on confirmation test sensitivity (mst 2, age 50, indolent_rate=0.2, confirmation test rate = 1)?
  dset_c_sens <- dset %>% filter(confirmation_test_rate==1)
  dset_c_sens <- dset_c_sens %>% select(-c(test_age,
                                           mean_sojourn_time,
                                           indolent_rate,
                                           confirmation_test_rate))
  plot_prospective_sensitivity(dset_c_sens,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)

  # compare different method of sensitivity
  plot_sens_compare(N=1000,
                    test_age=seq(40, 70, 2),
                    mean_sojourn_time=c(2, 5, 10),
                    ext=ext,
                    saveit=saveit)

  # compare different settings using analytic model
  dset <- expand_grid(test_age=c(50, 60, 70),
                      mean_sojourn_time,
                      indolent_rate=0,
                      confirmation_test_rate,
                      confirmation_test_sensitivity=1)
  plot_prospective_sensitivity_compare(N=N,
                                       dset=dset,
                                       start_age=start_age,
                                       preonset_rate=preonset_rate,
                                       onset_sensitivity=onset_sensitivity,
                                       clinical_sensitivity=clinical_sensitivity,
                                       significant_cancer_only=FALSE,
                                       ext=ext,
                                       saveit=saveit)
  
  # compare prospective sensitivity for progressive cancer only
  dset <- expand_grid(test_age=c(50, 60, 70),
                      mean_sojourn_time,
                      indolent_rate,
                      confirmation_test_rate=1,
                      confirmation_test_sensitivity=1)
  plot_prospective_sensitivity_compare(N=N,
                                       dset=dset,
                                       start_age=start_age,
                                       preonset_rate=preonset_rate,
                                       onset_sensitivity=onset_sensitivity,
                                       clinical_sensitivity=clinical_sensitivity,
                                       indolent_sensitivity=0.2,
                                       significant_cancer_only=TRUE,
                                       ext=ext,
                                       saveit=saveit)
  # compare prospective sensitivity under various settings
  dset <- expand_grid(test_age=c(50, 60, 70),
                      mean_sojourn_time,
                      indolent_rate,
                      confirmation_test_rate,
                      confirmation_test_sensitivity=1)
  plot_prospective_sensitivity_compare(N=N,
                                       dset=dset,
                                       start_age=start_age,
                                       preonset_rate=preonset_rate,
                                       onset_sensitivity=onset_sensitivity,
                                       clinical_sensitivity=clinical_sensitivity,
                                       indolent_sensitivity=0.2,
                                       significant_cancer_only=FALSE,
                                       ext=ext,
                                       saveit=saveit)
  
}
control(saveit=TRUE)

