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
library(rlang)
library(ggpubr)
library(ggpmisc)

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
# datestamp <- '2023-06-21'
# datestamp <- '2023-07-24'
# datestamp <- '2023-07-26'
# datestamp <- '2023-07-31'
# datestamp <- '2023-09-25'
# datestamp <- '2023-10-04'
# datestamp <- '2023-10-06'
# datestamp <- '2023-10-17'
# datestamp <- '2024-02-14'
#datestamp <- '2024-02-15'
#datestamp <- '2024-02-28'
#datestamp <- '2024-03-11'
#datestamp <- '2024-03-21'
#datestamp <- '2024-03-24'
datestamp <- '2024-03-25'

set.seed(12345)

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
                             length=20,
                             time=seq(0, sojourn_time, length=length),
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
    #tset <- tset %>%
    #  mutate(sensitivity=ifelse(sensitivity<=clinical_sensitivity, sensitivity, NA),
    #         sensitivity=ifelse(time<0|time>sojourn_time, NA, sensitivity))
  }else{
    tset <- tibble(time=time, 
                   sensitivity=onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/sojourn_time*time)))
  }
  return(tset)
}

plot_test_sensitivity <- function(sojourn_time,
                                  onset_sensitivity,
                                  clinical_sensitivity,
                                  ext='png',
                                  saveit=FALSE){
  sset0 <- tibble(sojourn_time)
  sset1 <- sset0 %>% mutate(sensitivity=map(sojourn_time,
                                            test_sensitivity,
                                            onset_sensitivity,
                                            clinical_sensitivity))
  sset1 <- sset1 %>% unnest(sensitivity)
  sset2 <- sset0 %>% mutate(sensitivity=map(sojourn_time,
                                            test_sensitivity,
                                            onset_sensitivity,
                                            clinical_sensitivity,
                                            is_indolent=TRUE))
  sset2 <- sset2 %>% unnest(sensitivity)
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               legend.position='bottom')
  gg <- ggplot(sset1)
  gg <- gg+geom_hline(aes(yintercept=clinical_sensitivity), linetype='dashed')
  # gg <- gg+geom_hline(aes(yintercept=onset_sensitivity), color='red')
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
  gg <- gg+geom_segment(data=sset1 %>% filter(sojourn_time == time),
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
  gg <- gg+geom_line(data=sset2,
                     aes(x=time,
                         y=sensitivity,
                         group=sojourn_time,
                         colour=sojourn_time),
                     linetype=2)
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

plot_test_sensitivity_overall_example <- function(onset_sensitivity=0,
                                                  clinical_sensitivity=0.8,
                                                  test_age=55){
  pset <- tibble(type='progressive',
                 onset_age=49,
                 sojourn_time=8,
                 clinical_age=onset_age+sojourn_time)
  pset <- pset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity,
                                          length=100))
  pset <- pset %>% unnest(sensitivity)
  iset <- tibble(type='indolent',
                 onset_age=49,
                 sojourn_time=8,
                 clinical_age=onset_age+sojourn_time)
  iset <- iset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity,
                                          is_indolent=TRUE,
                                          length=100))
  iset <- iset %>% unnest(sensitivity)
  #dset <- bind_rows(pset, iset)
  dset <- pset
  dset <- dset %>% mutate(age=onset_age+time)
  # identify prospective sensitivity within type
  xset <- dset %>% group_by(type)
  xset <- xset %>% filter(age == age[which.min(abs(age-test_age))])
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               plot.title=element_text(hjust=0.5),
               legend.position='none')
  gg <- ggplot(dset)
  gg <- gg+labs(title='Sensitivity across clinical stages')
  gg <- gg+geom_hline(yintercept=clinical_sensitivity,
                      color='#25848EFF',
                      linetype='dashed')
  gg <- gg+annotate('text',
                    x=100,
                    y=clinical_sensitivity,
                    label=sprintf('Diagnostic = %2.0f%%', 100*clinical_sensitivity),
                    color='#25848EFF',
                    hjust=1,
                    vjust=-0.5)
  gg <- gg+geom_line(aes(x=age,
                         y=sensitivity,
                         group=type,
                         linetype=type),
                     color='#25848EFF',
                     linewidth=0.35,
                     show.legend=TRUE)
  gg <- gg+geom_vline(xintercept=test_age,
                      linewidth=1.25,
                      color='orange',
                      alpha=0.75)
  gg <- gg+geom_segment(data=xset,
                        aes(x=60,
                            xend=age+1,
                            y=sensitivity,
                            yend=sensitivity),
                        color='#25848EFF',
                        arrow=arrow(length=unit(0.2, 'cm'), type='closed'),
                        show.legend=FALSE)
  gg <- gg+geom_text(data=xset,
                     aes(x=60,
                         y=sensitivity,
                         label=sprintf('Preclinical = %2.0f%%', 100*sensitivity)),
                     color='#25848EFF',
                     hjust=-0.1,
                     show.legend=FALSE)
  gg <- gg+scale_x_continuous(name='Age (years)',
                              limits=c(40, 100),
                              breaks=seq(40, 100, by=5),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='Sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.1),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_linetype_manual(name='', values=c(progressive='solid', indolent='dotted'))
  return(gg)
}

plot_test_sensitivity_stage_example <- function(onset_sensitivity=0,
                                                clinical_sensitivity=c(Early=0.3, Late=0.8),
                                                test_age=55){
  eset <- tibble(onset_age=48,
                 clinical_stage='Early',
                 sojourn_time=10,
                 clinical_age=onset_age+sojourn_time)
  eset <- eset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['Early'],
                                          length=100))
  eset <- eset %>% unnest(sensitivity)
  lset <- tibble(onset_age=52,
                 clinical_stage='Late',
                 sojourn_time=6,
                 clinical_age=onset_age+sojourn_time)
  lset <- lset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['Late'],
                                          length=100))
  lset <- lset %>% unnest(sensitivity)
  sset <- bind_rows(eset, lset)
  sset <- sset %>% mutate(age=onset_age+time)
  # identify prospective sensitivity within stage
  xset <- sset %>% group_by(clinical_stage)
  xset <- xset %>% filter(age == age[which.min(abs(age-test_age))])
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               plot.title=element_text(hjust=0.5),
               legend.position='none')
  gg <- ggplot(sset)
  gg <- gg+labs(title='Sensitivity within clinical stages')
  gg <- gg+geom_hline(data=xset,
                      aes(yintercept=clinical_sensitivity,
                          color=clinical_stage),
                      linetype='dashed')
  gg <- gg+geom_text(data=xset,
                     aes(x=100,
                         y=clinical_sensitivity,
                         color=clinical_stage,
                         label=sprintf('Diagnostic = %2.0f%%', 100*clinical_sensitivity)),
                     hjust=1,
                     vjust=-0.5,
                     show.legend=FALSE)
  gg <- gg+geom_line(aes(x=age,
                         y=sensitivity,
                         group=clinical_stage,
                         color=clinical_stage),
                     linewidth=0.35,
                     show.legend=TRUE)
  gg <- gg+geom_vline(xintercept=test_age,
                      linewidth=1.25,
                      color='orange',
                      alpha=0.75)
  gg <- gg+geom_segment(data=xset,
                        aes(x=60,
                            xend=age+1,
                            y=sensitivity,
                            yend=sensitivity,
                            color=clinical_stage),
                        arrow=arrow(length=unit(0.2, 'cm'), type='closed'),
                        show.legend=FALSE)
  gg <- gg+geom_text(data=xset,
                     aes(x=60,
                         y=sensitivity,
                         color=clinical_stage,
                         label=sprintf('Preclinical = %2.0f%%', 100*sensitivity)),
                     hjust=-0.1,
                     show.legend=FALSE)
  gg <- gg+scale_x_continuous(name='Age (years)',
                              limits=c(40, 100),
                              breaks=seq(40, 100, by=5),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='Sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.1),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_color_viridis(name='Clinical\nstage',
                               discrete=TRUE,
                               begin=0.65,
                               end=0.25)
  gg <- gg+guides(col=guide_legend(override.aes=list(linewidth=1, alpha=1)))
  return(gg)
}

plot_test_sensitivity_examples <- function(ext='png', saveit=FALSE){
  go <- plot_test_sensitivity_overall_example()
  gs <- plot_test_sensitivity_stage_example()
  gg <- ggarrange(go, gs, nrow=1)
  gg <- gg+bgcolor('white')
  print(gg)
  if(saveit){
    filename <- str_glue('prospective_sensitivity_examples_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=6,
           width=12)
  }
}

plot_test_sensitivity_stage <- function(sojourn_time=tibble(early=seq(2, 10, by=2),
                                                            late=seq(6, 14, by=2)),
                                        onset_sensitivity,
                                        clinical_sensitivity=c(early=0.4, late=0.8),
                                        ext='png',
                                        saveit=FALSE){
  eset <- sojourn_time %>% select(sojourn_time=early)
  eset <- eset %>% mutate(clinical_stage='early',
                          sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['early']))
  eset <- eset %>% unnest(sensitivity)
  lset <- sojourn_time %>% select(sojourn_time=late)
  lset <- lset %>% mutate(clinical_stage='late',
                          sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['late']))
  lset <- lset %>% unnest(sensitivity)
  sset <- bind_rows(eset, lset)
  cset <- sset %>% group_by(clinical_stage)
  cset <- cset %>% summarize(min=min(sojourn_time),
                             max=max(sojourn_time))
  cset <- cset %>% mutate(sensitivity=clinical_sensitivity[clinical_stage])
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               legend.position='bottom')
  gg <- ggplot(sset)
  gg <- gg+geom_hline(yintercept=clinical_sensitivity, linetype='dashed')
  gg <- gg+geom_line(aes(x=time,
                         y=sensitivity,
                         group=interaction(clinical_stage, sojourn_time),
                         color=clinical_stage),
                     alpha=0.9,
                     linewidth=0.75)
  gg <- gg+geom_segment(data=cset,
                        aes(x=min,
                            xend=max,
                            y=sensitivity,
                            yend=sensitivity),
                        alpha=0.9,
                        linewidth=2)
  gg <- gg+scale_x_continuous(name='Years since onset',
                              limits=c(0, max(sojourn_time)),
                              breaks=seq(0, max(sojourn_time), by=2),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='Sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.1),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_color_viridis(name='Clinical stage',
                               discrete=TRUE,
                               begin=0.65,
                               end=0.25)
  print(gg)
  if(saveit){
    filename <- str_glue('prospective_sensitivity_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=6,
           width=6)
  }
}

plot_test_sensitivity_stage_experiment <- function(mst=c(etoc=6, etol=4),
                                                   onset_sensitivity=0,
                                                   clinical_sensitivity=c(early=0.4, late=0.8),
                                                   size=1000,
                                                   ext='png',
                                                   saveit=FALSE){
  eset <- tibble(sojourn_time=rexp(size, rate=1/mst['etoc']))
  eset <- eset %>% mutate(clinical_stage='early',
                          sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['early']))
  eset <- eset %>% unnest(sensitivity)
  lset <- tibble(sojourn_time=rexp(size, rate=1/mst['etol']))
  lset <- lset %>% mutate(clinical_stage='late',
                          sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['late']))
  lset <- lset %>% unnest(sensitivity)
  sset <- bind_rows(eset, lset)
  cset <- sset %>% group_by(clinical_stage)
  cset <- cset %>% summarize(min=min(sojourn_time),
                             max=max(sojourn_time))
  cset <- cset %>% mutate(sensitivity=clinical_sensitivity[clinical_stage])
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               legend.position='bottom')
  gg <- ggplot(sset)
  gg <- gg+geom_hline(yintercept=clinical_sensitivity, linetype='dashed')
  gg <- gg+geom_line(aes(x=time,
                         y=sensitivity,
                         group=interaction(clinical_stage, sojourn_time),
                         color=clinical_stage),
                     alpha=0.5,
                     linewidth=0.25)
  gg <- gg+geom_segment(data=cset,
                        aes(x=min,
                            xend=max,
                            y=sensitivity,
                            yend=sensitivity),
                        alpha=0.9,
                        linewidth=1.5)
  gg <- gg+scale_x_continuous(name='Years since onset',
                              limits=c(0, max(cset$max)),
                              breaks=seq(0, max(cset$max), by=2))
  gg <- gg+scale_y_continuous(name='Sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.1),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_color_viridis(name='Clinical stage',
                               discrete=TRUE,
                               begin=0.65,
                               end=0.25)
  gg <- gg+guides(col=guide_legend(override.aes=list(linewidth=1, alpha=1)))
  print(gg)
  if(saveit){
    filename <- str_glue('prospective_sensitivity_experiment_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=6,
           width=6)
  }
}

plot_test_sensitivity_age_stage_experiment <- function(dset,
                                                       onset_sensitivity=0,
                                                       clinical_sensitivity=c(Early=0.3, Late=0.8),
                                                       test_age=55,
                                                       age=40,
                                                       mht=10,
                                                       size=400,
                                                       ext='png',
                                                       saveit=FALSE){
  oset <- tibble(id=seq(size),
                 onset_age=age+rexp(size, rate=1/mht),
                 clinical_stage=rbinom(n=size,
                                       size=1,
                                       prob=with(dset, Proportion[Stage == 'Early'])))
  oset <- oset %>% mutate(clinical_stage=factor(clinical_stage,
                                                levels=c(1, 0),
                                                labels=c('Early', 'Late')))
  eset <- oset %>% filter(clinical_stage == 'Early')
  eset <- eset %>% mutate(sojourn_time=rexp(length(id), rate=1/with(dset, MST[Stage == 'Early'])),
                          clinical_age=onset_age+sojourn_time)
  eset <- eset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['Early']))
  eset <- eset %>% unnest(sensitivity)
  lset <- oset %>% filter(clinical_stage == 'Late')
  lset <- lset %>% mutate(sojourn_time=rexp(length(id), rate=1/with(dset, MST[Stage == 'Late'])),
                          clinical_age=onset_age+sojourn_time)
  lset <- lset %>% mutate(sensitivity=map(sojourn_time,
                                          test_sensitivity,
                                          onset_sensitivity,
                                          clinical_sensitivity['Late']))
  lset <- lset %>% unnest(sensitivity)
  sset <- bind_rows(eset, lset)
  sset <- sset %>% mutate(age=onset_age+time)
  sset <- sset %>% group_by(id)
  sset <- sset %>% mutate(detected=onset_age <= test_age & test_age < clinical_age)
  # calculate empirical sensitivity within stage
  xset <- sset %>% filter(detected)
  xset <- xset %>% filter(age == age[which.min(abs(age-test_age))])
  xset <- xset %>% ungroup()
  xset <- xset %>% group_by(Stage=clinical_stage)
  xset <- xset %>% summarize(Sensitivity=mean(sensitivity))
  tset <- full_join(dset, xset, by='Stage')
  tset <- tset %>% mutate(Proportion=sprintf('%4.2f', Proportion),
                          Sensitivity=sprintf('%4.2f', Sensitivity))
  gt_theme <- gridExtra::ttheme_minimal()
  cols <- c('#2FB47CFF', '#3B528BFF')
  gt_theme$core$fg_params$col <- cols
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'),
               legend.position='bottom')
  gg <- ggplot(sset)
  gg <- gg+geom_hline(yintercept=clinical_sensitivity, linetype='dashed')
  gg <- gg+geom_line(aes(x=age,
                         y=sensitivity,
                         group=interaction(id, clinical_stage),
                         color=clinical_stage,
                         alpha=detected),
                     linewidth=0.35,
                     show.legend=c(color=TRUE, alpha=FALSE))
  gg <- gg+geom_vline(xintercept=test_age,
                      linewidth=1.25,
                      color='orange',
                      alpha=0.75)
  gg <- gg+annotate(geom='table',
                    x=100,
                    y=1,
                    vjust=1,
                    hjust=1,
                    label=list(tset),
                    table.theme=gt_theme)
  gg <- gg+scale_x_continuous(name='Age (years)',
                              limits=c(40, 100),
                              breaks=seq(40, 100, by=5),
                              expand=c(0, 0))
  gg <- gg+scale_y_continuous(name='Sensitivity',
                              limits=c(0, 1),
                              breaks=seq(0, 1, by=0.1),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  gg <- gg+scale_color_viridis(name='Clinical stage',
                               discrete=TRUE,
                               begin=0.65,
                               end=0.25)
  gg <- gg+scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.2))
  gg <- gg+guides(col=guide_legend(override.aes=list(linewidth=1, alpha=1)))
  print(gg)
  if(saveit){
    filename <- str_glue('prospective_sensitivity_age_stage_{svec["etol"]}_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=5,
           width=8)
  }
  return(gg)
}

plot_test_sensitivity_age_stage_experiment_table <- function(ext='png', saveit=FALSE){
  mset <- expand_grid(Early=c(2, 10),
                      Late=c(2, 10))
  mset <- mset %>% mutate(Scenario=seq(n()))
  mset <- mset %>% pivot_longer(cols=-Scenario,
                                names_to='Stage',
                                values_to='MST')
  mset <- mset %>% mutate(Proportion=case_when(Scenario == 1 & Stage == 'Early' ~ 0.50,
                                               Scenario == 2 & Stage == 'Early' ~ 0.83,
                                               Scenario == 3 & Stage == 'Early' ~ 0.17,
                                               Scenario == 4 & Stage == 'Early' ~ 0.50))
  mset <- mset %>% group_by(Scenario)
  mset <- mset %>% mutate(Proportion=if_else(Stage == 'Early',
                                             Proportion,
                                             1-sum(Proportion, na.rm=TRUE)))
  mset <- mset %>% nest(data=-Scenario)
  mset <- mset %>% mutate(plot=map(data, plot_test_sensitivity_age_stage_experiment))
  mset <- mset %>% mutate(Scenario=factor(Scenario, levels=c(2, 4, 1, 3)))
  mset <- mset %>% arrange(Scenario)
  gg <- ggarrange(plotlist=mset$plot,
                  nrow=2,
                  ncol=2,
                  common.legend=TRUE,
                  legend='bottom')
  gg <- gg+bgcolor('white')
  print(gg)
  if(saveit){
    filename <- str_glue('prospective_sensitivity_age_stage_experiments_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           height=10,
           width=12)
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
                                         confirmation_test_rate=NULL,
                                         confirmation_test_sensitivity=NULL){
  # simulate ages at preclinical onset and clinical diagnosis
  dset <- tibble(onset_age=start_age+rexp(N, rate=preonset_rate),
                 sojourn_time=rexp(N, rate=1/mean_sojourn_time),
                 indolent=rbinom(N, 1, indolent_rate)==1,
                 clinical_age=ifelse(indolent, Inf, onset_age+sojourn_time))
  # confirm retrospective test sensitivity at clinical diagnosis
  cset <- dset %>% filter(is.finite(clinical_age))
  cset <- cset %>% group_by(sojourn_time, indolent) %>%
    mutate(retro_sens=test_sensitivity(sojourn_time=sojourn_time,
                                                      onset_sensitivity,
                                                      clinical_sensitivity,
                                                      is_indolent=indolent,
                                                      time=sojourn_time)$sensitivity) %>%
    ungroup()
  cset <- cset %>% mutate(retro_sens=ifelse(is.na(retro_sens), clinical_sensitivity, retro_sens))
  # ccheck <- cset %>% with(all(retro_sens == clinical_sensitivity))
  # ccheck %>% stopifnot()
  # calculate true prospective test sensitivity among preclinical at test age
  pset <- dset %>% filter(onset_age <= !!test_age & !!test_age < clinical_age)
  pset <- pset %>% mutate(prosp_sens_all=test_sensitivity(sojourn_time=sojourn_time,
                                                      onset_sensitivity,
                                                      clinical_sensitivity,
                                                      time=test_age-onset_age)$sensitivity)
  # update indolent sensitivity
  pset <- pset %>% mutate(prosp_sens_all=case_when(indolent==1~test_sensitivity(sojourn_time=sojourn_time,
                                                                                onset_sensitivity=onset_sensitivity,
                                                                                clinical_sensitivity=clinical_sensitivity,
                                                                                is_indolent=TRUE,
                                                                                time=(test_age-onset_age))$sensitivity,
                                                   TRUE~prosp_sens_all))
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
## f(u, preonset_rate): function of preclinical onset time and onset rate
##            u: time interval between initial age to preclinical onset
## g(s, mean_sojourn_time): function of sojourn time
##            s: sojourn time, time between preclinical onset and clinical onset
##            mean_sojourn_time: mean sojourn time
## G(t, u, t0, mean_sojourn_time): cumulative distribution function of sojourn time and mean sojourn time
##            t: screen time
##            u: time interval between initial age to preclinical onset
##            t0: initial time
##            mean_sojourn_time: mean sojourn time
## h(t, s, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): true sensitivity for both cases
##            t: screen time
##            s: sojourn time, time between preclinical onset and clinical onset
##            u: time interval between initial age to preclinical onset
##            t0: initial time
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## h_g(t, s, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): function of true sensitivity and sojourn time for progressive cases
##            t: screen time
##            u: preclinical onset age
##            s: clinical onset age
##            t0: initial time
##            onset_sensitivity: sensitivity at prelinical onset, lower bound of sensitivity
##            clinical_sensitivity: sensitivity at clinical onset, upper bound of sensitivity
##            mean_sojourn_time: mean sojourn time
## integral_h_g_progressive(u, t, s, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time): integral of f_g for progressive cases
##            u: preclinical onset age
##            t: screen time
##            s: clinical onset age
##            t0: initial time
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

f <- function(u, preonset_rate){
  preonset_rate*exp(-preonset_rate*u)
}

g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

G <- function(t, u, t0, mean_sojourn_time){
  tp <- t0+u # preclinical onset time
  1-exp(-(1/mean_sojourn_time)*(t-tp))
}

# sensitivity function
h <- function(t, s, u, t0, onset_sensitivity, clinical_sensitivity){
  tp <- t0+u # preclinical onset time
  # progressive case
  alpha <- onset_sensitivity
  beta <- (clinical_sensitivity-onset_sensitivity)/s
  sensitivity <- alpha+beta*(t-tp)
  sensitivity_progressive <- ifelse(sensitivity<=clinical_sensitivity, sensitivity, NA)
  # indolent case
  sensitivity_indolent <- onset_sensitivity+(clinical_sensitivity-onset_sensitivity)*(1-exp(-1/s*(t-tp)))
  
  return(c(sensitivity_progressive, sensitivity_indolent))
}

f_h <- function(u, s, t0, t, preonset_rate, onset_sensitivity, clinical_sensitivity){
  f(u, preonset_rate)*h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)[2]
}
  
h_g <- function(s, t, u, t0, onset_sensitivity, clinical_sensitivity, mean_sojourn_time, indolent=FALSE){
  if(indolent){
    return(h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)[2]*g(s, mean_sojourn_time))
  }
  return(h(t, s, u, t0, onset_sensitivity, clinical_sensitivity)[1]*g(s, mean_sojourn_time))
}

f_G <- function(u, t, t0, preonset_rate, mean_sojourn_time){
  f(u, preonset_rate)*(1-G(t, u, t0, mean_sojourn_time))
}

integral_h_g <- function(u, t, t0, preonset_rate, onset_sensitivity, clinical_sensitivity, mean_sojourn_time, indolent){
  tp <- t0+u # preclinical onset time
  if(indolent){
    f(u, preonset_rate)*
      integrate(Vectorize(h_g),
                lower=0, upper=Inf,
                t=t,
                u=u,
                t0=t0,
                onset_sensitivity=onset_sensitivity,
                clinical_sensitivity=clinical_sensitivity,
                mean_sojourn_time=mean_sojourn_time,
                indolent=indolent)$value
  }else{
    f(u, preonset_rate)*
      integrate(Vectorize(h_g),
                lower=t-tp, upper=Inf,
                t=t,
                u=u,
                t0=t0,
                onset_sensitivity=onset_sensitivity,
                clinical_sensitivity=clinical_sensitivity,
                mean_sojourn_time=mean_sojourn_time,
                indolent=indolent)$value
  }
}

prospective_sens_num <- function(t, t0,
                                 preonset_rate,
                                 onset_sensitivity,
                                 clinical_sensitivity,
                                 mean_sojourn_time,
                                 indolent_rate){
  prospective_sens_num_integral_progressive <- (1-indolent_rate) *
    integrate(Vectorize(integral_h_g),
              lower=0, upper=t-t0,
              t=t,
              t0=t0,
              preonset_rate=preonset_rate,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time,
              indolent=FALSE)$value
  prospective_sens_num_integral_indolent <- indolent_rate *
    integrate(Vectorize(integral_h_g),
              lower=0, upper=t-t0,
              t=t,
              t0=t0,
              preonset_rate=preonset_rate,
              onset_sensitivity=onset_sensitivity,
              clinical_sensitivity=clinical_sensitivity,
              mean_sojourn_time=mean_sojourn_time,
              indolent=TRUE)$value
  
  prospective_sens_num_integral_progressive+prospective_sens_num_integral_indolent
}

prospective_sens_denom <- function(t, t0, preonset_rate, mean_sojourn_time, indolent_rate){
  prospective_sens_denom_integral_progressive <- (1-indolent_rate)*
    integrate(Vectorize(f_G),
              lower=0, upper=t-t0,
              t=t, 
              t0=t0,
              preonset_rate=preonset_rate, 
              mean_sojourn_time=mean_sojourn_time)$value
  prospective_sens_denom_integral_indolent <- indolent_rate*
    integrate(Vectorize(f),
              lower=0, upper=t-t0,
              preonset_rate=preonset_rate)$value
  
  prospective_sens_denom_integral_progressive+prospective_sens_denom_integral_indolent
}


prospective_sens_analyt <- function(start_age,
                             test_age,
                             preonset_rate,
                             mean_sojourn_time,
                             onset_sensitivity,
                             clinical_sensitivity,
                             indolent_rate,
                             confirmation_test_rate=1,
                             confirmation_test_sensitivity=1){
  
  numerator <- prospective_sens_num(t=test_age, 
                                    t0=start_age,
                                    preonset_rate=preonset_rate,
                                    onset_sensitivity=onset_sensitivity, 
                                    clinical_sensitivity=clinical_sensitivity, 
                                    mean_sojourn_time=mean_sojourn_time,
                                    indolent_rate=indolent_rate)
  
  denomerator <- prospective_sens_denom(t=test_age, 
                                        t0=start_age, 
                                        preonset_rate=preonset_rate, 
                                        mean_sojourn_time=mean_sojourn_time, 
                                        indolent_rate=indolent_rate)    
    
  return(confirmation_test_rate*confirmation_test_sensitivity*numerator/denomerator)
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
                         start_age,
                         test_age,
                         preonset_rate,
                         mean_sojourn_time,
                         onset_sensitivity,
                         clinical_sensitivity,
                         indolent_rate,
                         confirmation_test_rate,
                         confirmation_test_sensitivity
  ){
  analytic_sens <- prospective_sens_analyt(start_age=start_age,
                                    test_age=test_age,
                                    preonset_rate=preonset_rate,
                                    mean_sojourn_time=mean_sojourn_time,
                                    onset_sensitivity=onset_sensitivity,
                                    clinical_sensitivity=clinical_sensitivity,
                                    indolent_rate=indolent_rate,
                                    confirmation_test_rate=confirmation_test_rate,
                                    confirmation_test_sensitivity=confirmation_test_sensitivity)
  
  sim_sens_tb <- prospective_test_sensitivity(N=N,
                                              start_age=start_age,
                                              preonset_rate=preonset_rate,
                                              mean_sojourn_time=mean_sojourn_time,
                                              indolent_rate=indolent_rate,
                                              test_age=test_age,
                                              onset_sensitivity=onset_sensitivity,
                                              clinical_sensitivity=clinical_sensitivity,
                                              confirmation_test_rate=confirmation_test_rate,
                                              confirmation_test_sensitivity=confirmation_test_sensitivity)[[1]]
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
plot_sens_compare <- function(N=1000,
                              start_age=40,
                              test_age=seq(40, 70, 2),
                              preonset_rate=0.1,
                              mean_sojourn_time=c(2, 5, 10),
                              onset_sensitivity=0.2,
                              clinical_sensitivity=0.8,
                              indolent_rate=c(0, 0.2, 0.6),
                              confirmation_test_rate=c(0.6, 0.8, 1),
                              confirmation_test_sensitivity=c(0.4, 0.6, 1),
                              # N,
                              # start_age,
                              # test_age,
                              # preonset_rate,
                              # mean_sojourn_time,
                              # onset_sensitivity,
                              # clinical_sensitivity,
                              # indolent_rate,
                              # confirmation_test_rate=1,
                              # confirmation_test_sensitivity=1,
                              ext='png',
                              saveit=FALSE){
  dset <- expand_grid(test_age, mean_sojourn_time, indolent_rate, confirmation_test_rate, confirmation_test_sensitivity)
  dset <- dset %>% group_by(test_age, mean_sojourn_time, indolent_rate, confirmation_test_rate, confirmation_test_sensitivity)
  dset <- dset %>% mutate(results=sens_compare(N=N,
                                               start_age=start_age,
                                               test_age=test_age,
                                               preonset_rate=preonset_rate,
                                               mean_sojourn_time=mean_sojourn_time,
                                               onset_sensitivity=onset_sensitivity,
                                               clinical_sensitivity=clinical_sensitivity,
                                               indolent_rate=indolent_rate,
                                               confirmation_test_rate=confirmation_test_rate,
                                               confirmation_test_sensitivity=confirmation_test_sensitivity))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  dset <- dset %>% pivot_longer(cols=ends_with('sens'),
                                names_to='type', values_to='sensitivity')
  # dset$test_age <- factor(dset$test_age)
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  # compare indolent rate and mean sojourn time
  gg <- dset %>% filter(confirmation_test_rate==1, confirmation_test_sensitivity==1) %>% 
    ggplot(aes(x=test_age, y=sensitivity, color=type))
  gg <- gg+geom_line()
  gg <- gg+facet_grid(indolent_rate~mean_sojourn_time, labeller=label_both)
  gg <- gg+xlab('Test age')
  gg <- gg+ylab('Sensitivity')
  gg <- gg+ylim(0,1)
  gg <- gg+guides(color=guide_legend(title='Types of prospective sensitivity'))
  gg <- gg+scale_color_discrete(labels=c('Analytic', 'Simulation'))
  gg <- gg+theme(legend.position='bottom')
  print(gg)
  
  if((length(confirmation_test_rate)+length(confirmation_test_sensitivity))>2){
    gg2 <- dset %>% filter(indolent_rate==0, mean_sojourn_time==2) %>%
      ggplot(aes(x=test_age, y=sensitivity, color=type))
    # gg <- gg+geom_boxplot()
    gg2 <- gg2+geom_line()
    gg2 <- gg2+facet_grid(confirmation_test_rate~confirmation_test_sensitivity, labeller=label_both)
    gg2 <- gg2+xlab('Test age')
    gg2 <- gg2+ylab('Sensitivity')
    gg2 <- gg2+ylim(0,1)
    gg2 <- gg2+guides(color=guide_legend(title='Types of prospective sensitivity'))
    gg2 <- gg2+scale_color_discrete(labels=c('Analytic', 'Simulation'))
    gg2 <- gg2+theme(legend.position='bottom')
    print(gg2)
  }
  # compare confirmation test

  if(saveit){
    filename <- str_glue('fig_sens_compare_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           width=6,
           height=5)
    if((length(confirmation_test_rate)+length(confirmation_test_sensitivity))>2){
      filename2 <- str_glue('fig_sens_confirmation_compare_{datestamp}.{ext}')
      ggsave(here('plot', filename2),
             plot=gg2,
             width=6,
             height=5)
    }
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
                                                 method='analytic',
                                                 ext='png',
                                                 saveit=FALSE){
  # analytic_sens <- prospective_sens_analyt(start_age, test_age, onset_sensitivity, clinical_sensitivity, preonset_rate, mean_sojourn_time)
  dset <- dset %>% group_by(across(all_of(names(dset))))
  if(method=='simulation'){
    dset <- dset %>% mutate(prosp_sens_all=prospective_test_sensitivity(N=N,
                                                                        start_age=start_age,
                                                                        preonset_rate=preonset_rate,
                                                                        mean_sojourn_time=mean_sojourn_time,
                                                                        indolent_rate=indolent_rate,
                                                                        test_age=test_age,
                                                                        onset_sensitivity=onset_sensitivity,
                                                                        clinical_sensitivity=clinical_sensitivity,
                                                                        confirmation_test_rate=confirmation_test_rate,
                                                                        confirmation_test_sensitivity=confirmation_test_sensitivity))
    dset <- dset %>% mutate(mean_prosp_sens=map(prosp_sens_all, ~mean(.x$prosp_sens_all)))
  }
  if(method=='analytic'){
    dset <- dset %>% mutate(prosp_sens_all=
                              prospective_sens_analyt(start_age=start_age,
                                                      test_age=test_age,
                                                      preonset_rate=preonset_rate,
                                                      mean_sojourn_time=mean_sojourn_time,
                                                      onset_sensitivity=onset_sensitivity,
                                                      clinical_sensitivity=clinical_sensitivity,
                                                      indolent_rate=indolent_rate,
                                                      confirmation_test_rate=confirmation_test_rate,
                                                      confirmation_test_sensitivity=confirmation_test_sensitivity)
    )
    dset <- dset %>% mutate(mean_prosp_sens=prosp_sens_all)
  }
  ylabel <- 'Prospective sensitivity for any cancer'
  
  dset <- dset %>% select(-prosp_sens_all)
  dset <- dset %>% ungroup()
  dset <- dset %>% unnest(mean_prosp_sens)
  dset <- dset %>% mutate(degradation=clinical_sensitivity-mean_prosp_sens,
                          test_age=factor(test_age),
                          mean_sojourn_time=factor(mean_sojourn_time),
                          indolent_rate=paste('Indolent\nfraction', sprintf('%2.0f%%', 100*indolent_rate)),
                          indolent_rate=factor(indolent_rate)
                          # confirmation_test_sensitivity=paste('Confirmation test\nsensitivity', sprintf('%2.0f%%', 100*confirmation_test_sensitivity)),
                          # confirmation_test_sensitivity=factor(confirmation_test_sensitivity,
                          #                                     levels=c('Confirmation test\nsensitivity 60%',
                          #                                              'Confirmation test\nsensitivity 80%',
                          #                                              'Confirmation test\nsensitivity 100%'))
                          )
  dset <- dset %>% mutate( 
    confirmation_test_rate=factor(paste('Confirmation test\nfrequency/sensitivity', 
                                        sprintf('%2.0f%%', 100*confirmation_test_rate)),
                                  levels=paste('Confirmation test\nfrequency/sensitivity', 
                                               sprintf('%2.0f%%', 100*sort(unique(dset$confirmation_test_rate), decreasing=FALSE)))
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
                    position=position_dodge2(width=1))
  gg <- gg+geom_linerange(aes(#xmin=test_age,
                              ymin=mean_prosp_sens,
                              x=test_age,
                              #y=mean_prosp_sens,
                              #xmax=test_age,
                              ymax=clinical_sensitivity,
                              color=mean_sojourn_time),
                          alpha=0.3,
                          position=position_dodge2(width=1))
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
  gg <- gg+facet_grid(indolent_rate~confirmation_test_rate)
  gg <- gg+scale_x_discrete(name='Age at screening test (years)',
                            expand=c(0, 0))
  gg <- gg+scale_y_continuous(limits=c(0, 1),
                              breaks=seq(0, 1, by=0.2),
                              labels=label_percent(accuracy=1),
                              expand=c(0, 0))
  # gg <- gg+scale_color_viridis(name='Mean sojourn\ntime (years)', discrete=TRUE)
  # gg <- gg+scale_fill_viridis(name='Mean sojourn\ntime (years)', discrete=TRUE)
  gg <- gg+scale_color_manual(name = 'Mean sojourn\ntime (years)',
                             values = setNames(c("grey10", "grey40", "grey70"), 
                                               levels(dset$mean_sojourn_time)))
  gg <- gg+scale_fill_manual(name = 'Mean sojourn\ntime (years)',
                             values = setNames(c("grey10", "grey40", "grey70"), 
                                               levels(dset$mean_sojourn_time)))
  print(gg)
  if(saveit){
    filename <- str_glue('plot_compare_setting_{method}_{datestamp}.{ext}')
    ggsave(here('plot', filename),
           plot=gg,
           width=10,
           height=5)
  }
}

##################################################
# Estimate confidence interval of prospective sensitivity estimates from simulation
##################################################
summary_stats <- function(dset, 
                         start_age, 
                         onset_sensitivity,
                         clinical_sensitivity,
                         preonset_rate,
                         alpha){
  # summary statistics
  sset <- 
    dset %>%
    group_by(test_age, mean_sojourn_time, indolent_rate) %>%
    summarise(mean_sens=mean(prosp_sens_all),
              sd_sens=sd(prosp_sens_all),
              n=n())
  # estimate standard error
  sset <- sset %>% mutate(se=sd_sens/sqrt(n))
  # find the critical value for 95% confidence level
  sset <- sset %>% mutate(critical=ifelse(n<30, qt(1-alpha/2, df=n-1), qnorm(1-alpha/2)))
  # calculate the margin of error
  sset <- sset %>% mutate(margin_of_error=critical*se)
  # calculate the confidence interval
  sset <- sset %>% mutate(ci_lower=mean_sens-margin_of_error)
  sset <- sset %>% mutate(ci_upper=mean_sens+margin_of_error)
  
  # Calculate the confidence interval
  sset <- sset %>% 
    group_by(test_age, mean_sojourn_time, indolent_rate) %>% 
    mutate(prosp_sens=prospective_sens_analyt(start_age=!!start_age,
                                test_age=test_age,
                                onset_sensitivity=!!onset_sensitivity,
                                clinical_sensitivity=!!clinical_sensitivity,
                                preonset_rate=!!preonset_rate,
                                mean_sojourn_time=mean_sojourn_time,
                                indolent_rate=indolent_rate))
  return(sset)
}

summary_plot <- function(dset, 
                         start_age, 
                         onset_sensitivity,
                         clinical_sensitivity,
                         preonset_rate,
                         alpha=0.05,
                         ext,
                         saveit){
  sset <- summary_stats(dset,
                start_age, 
                onset_sensitivity,
                clinical_sensitivity,
                preonset_rate,
                alpha)
  
  label_test_age <- paste0("Test age: ", unique(sset$test_age))
  names(label_test_age) <- unique(sset$test_age)
  label_indolent_rate <- paste0("Indolent rate: ", percent(unique(sset$indolent_rate)))
  names(label_indolent_rate) <- unique(sset$indolent_rate)
  gg1 <- 
    ggplot(sset, aes(x=factor(mean_sojourn_time), y=prosp_sens, color=n)) +
    geom_point() +  # Plot the mean values
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0.2) +  # Add error bars
    facet_grid(test_age~indolent_rate, labeller=labeller(test_age=label_test_age,
                                                         indolent_rate=label_indolent_rate)) +
    scale_color_gradient(name="Group size", 
                         low="lightgrey", 
                         high="black") +
    labs(y="Prospective sensitivity", x="Mean sojourn time") +
    theme_update(strip.background=element_rect(colour=NA, fill=NA),
                 strip.text=element_text(size=10)) +
    theme_set(theme_classic())
  print(gg1)
  
  if(saveit){
    filename <- str_glue('fig_est_ci_{datestamp}.{ext}')
    ggsave(plot=print(gg1), here('plot', filename), width=6, height=5)
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
                    ext='png',
                    saveit=FALSE){
  # visualize how test sensitivity increases over sojourn time
  #plot_test_sensitivity(sojourn_time=seq(2, 10, by=2),
  #                      onset_sensitivity,
  #                      clinical_sensitivity,
  #                      ext=ext,
  #                      saveit=saveit)

  # visualize how test sensitivity depends on clinical stage
  #plot_test_sensitivity_stage(sojourn_time=tibble(early=seq(2, 6),
  #                                                late=seq(6, 10)),
  #                            onset_sensitivity=0,
  #                            clinical_sensitivity=c(early=0.3, late=0.8),
  #                            ext=ext,
  #                            saveit=saveit)

  # visualize how test sensitivity depends on clinical stage
  #plot_test_sensitivity_stage_experiment(mst=c(etoc=6, etol=4),
  #                                       onset_sensitivity=0,
  #                                       clinical_sensitivity=c(early=0.3, late=0.8),
  #                                       ext=ext,
  #                                       saveit=saveit)

  # visualize conceptual model of sensitivity
  plot_test_sensitivity_examples(ext=ext, saveit=saveit)

  # visualize how test sensitivity depends on age and clinical stage
  plot_test_sensitivity_age_stage_experiment_table(ext=ext, saveit=saveit)

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
                                                               confirmation_test_rate=confirmation_test_rate,
                                                               confirmation_test_sensitivity=confirmation_test_sensitivity)[[1]])
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
  # how does bias depend on confirmation test (mst 2, age 50, indolent_rate=0.2, confirmation test sensitivity=1)?
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
                                                               confirmation_test_rate,
                                                               confirmation_test_sensitivity))
  dset <- dset %>% unnest(results)
  dset <- dset %>% ungroup()
  ## how does bias depend on confirmation test rate (mst 2, age 50, indolent_rate=0.2, confirmation test sensitivity=1)?
  dset_c_rate <- dset %>% filter(confirmation_test_sensitivity==1)
  dset_c_rate <- dset_c_rate %>% select(-c(test_age,
                                           mean_sojourn_time,
                                           indolent_rate,
                                           confirmation_test_sensitivity))
  plot_prospective_sensitivity(dset_c_rate,
                               clinical_sensitivity,
                               ext=ext,
                               saveit=saveit)
  ## how does bias depend on confirmation test sensitivity (mst 2, age 50, indolent_rate=0.2, confirmation test rate=1)?
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
                    start_age=40,
                    test_age=seq(40, 70, 2),
                    preonset_rate=0.1,
                    mean_sojourn_time=c(2, 5, 10),
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    indolent_rate=c(0, 0.2, 0.6),
                    confirmation_test_rate=1,
                    confirmation_test_sensitivity=1,
                    ext=ext,
                    saveit=saveit)
  ## adding confirmation test
  plot_sens_compare(N=1000,
                    start_age=40,
                    test_age=seq(40, 70, 2),
                    preonset_rate=0.1,
                    mean_sojourn_time=c(2, 5, 10),
                    onset_sensitivity=0.2,
                    clinical_sensitivity=0.8,
                    indolent_rate=c(0, 0.2, 0.6),
                    confirmation_test_rate=c(0.6, 0.8, 1),
                    confirmation_test_sensitivity=c(0.4, 0.6, 1),
                    ext=ext,
                    saveit=saveit)
  
  # compare different scenarios using analytic or simulation method
  dset <- expand_grid(test_age=c(50, 60, 70),
                      mean_sojourn_time=c(2, 5, 10),
                      indolent_rate=c(0, 0.2, 0.6),
                      confirmation_test_rate=c(0.6, 0.8, 1),
                      confirmation_test_sensitivity=1)
  ## analytic method
  plot_prospective_sensitivity_compare(N=N,
                                       dset=dset,
                                       start_age=40,
                                       preonset_rate=0.1,
                                       onset_sensitivity=0.2,
                                       clinical_sensitivity=0.8,
                                       method='analytic',
                                       ext=ext,
                                       saveit=saveit)
  ## simulation
  plot_prospective_sensitivity_compare(N=1000,
                                       dset=dset,
                                       start_age=40,
                                       preonset_rate=0.1,
                                       onset_sensitivity=0.2,
                                       clinical_sensitivity=0.8,
                                       method='simulation',
                                       ext=ext,
                                       saveit=saveit)
  # Plot estimates and 95% CI
  summary_plot(dset, 
               start_age, 
               onset_sensitivity,
               clinical_sensitivity,
               preonset_rate,
               alpha=0.05,
               ext,
               saveit)
}
control(saveit=TRUE)

