##################################################
# Compare true sensitivity VS empirical sensitivity in the simulation
##################################################
library(ggplot2)
source("simulation_functions_empirical_sensitivity_expanded.R")
# datestamp <- '2023-10-16'
datestamp <- '2023-10-17'

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
      mutate(sensitivity=ifelse(sensitivity<=clinical_sensitivity, sensitivity, clinical_sensitivity))
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
rateexample=matrix(c(-.002,.002,0,0,-.16,.16,0,0,0),byrow=T,nrow=3)
#simulation method
sens_sim <- empirical.sensitivity.simulation(rate.matrix=rateexample,
                                start.dist=c(1,0,0),
                                screen.times=1:10,
                                post.screen.lookout=1,
                                clinical.cancer.state=3,
                                pre.clinical.cancer.state=2,
                                nreps=10000)
sens_sim2$type <- "Empirical sensitivity"
sens_sim2 <- sens_sim %>% pivot_longer(c("empirical_sensitivity", "mean_true_sensitivity"),
                                       names_to = "type",
                                       values_to = "sensitivity")
sens_sim2 <- sens_sim %>% select(time=obs.times, sensitivity=empirical_sensitivity)
sens_true <- test_sensitivity(sojourn_time=1/rateexample[2,3],
                              onset_sensitivity=0.2,
                              clinical_sensitivity=0.8,
                              is_indolent=FALSE,
                              time=1:10,
                              method='linear')
sens_true$type <- "True sensitivity"
sens_all <- rbind(sens_sim2, sens_true)

sens_all %>%
  ggplot(aes(x=time, y=sensitivity, color=type)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  ylim(0, 1)
ggsave(here("plot", "fig_ES_TS.png"))
