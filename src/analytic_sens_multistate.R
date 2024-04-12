library(here)
library(foreach)
library(doParallel)
library(ggplot2)
library(scales)
library(dplyr)
library(stringr)

source(here("src", "analytic_sens_multistate_funs.R"))

# datastamp <- '2024-02-26'
# datastamp <- '2024-02-28'
# datastamp <- '2024-03-15'
# datastamp <- '2024-03-21'
datastamp <- '2024-03-29'

# plot prospective sensitivity for test @ between early stage and clinical onset
fig_prosp_sens <- function(test_time, rate, m1, m2, m3 = NULL,
                           trans_type,
                           ext, saveit){
  ############
  # test_time: test time at t
  # rate: preclinical onset rate
  # m1: mean sojourn time for transition from early stage to clinical onset
  # m2: mean sojourn time for transition from early stage to late stage
  # m3: mean sojourn time for transition from late stage to clinical onset
  # trans_type: 1. pre to clinical only; 2. Pre to clinical or late stage; 3. 1&2 + test at late stage
  ############
  no_cores <- detectCores() - 1  # leaving one core free
  registerDoParallel(cores=no_cores)
  cset <- expand.grid(test_time = test_time, m1 = m1, m2 = m2)
  cset$i <- seq(nrow(cset))
  
  theme_set(theme_classic())
  theme_update(axis.ticks.length=unit(0.2, 'cm'))
  
  # Case 1: test positive among early stage to clinical
  if(trans_type == 1){
    prosp_sens_es_cl <- foreach(i = cset$i, .combine = "rbind") %dopar% {
      test_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i]) / 
        case_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i])
    }
    dset_es_cl <- data.frame(cset, prosp_sens = c(prosp_sens_es_cl))
    # initial age is 40
    dset_es_cl <- dset_es_cl %>% mutate(test_time = (40 + test_time))
    
    # ggplot(data = dset_es_cl, aes(x = test_time, y = prosp_sens, color = factor(m1))) +
    #   geom_line() +
    #   facet_grid(. ~ m2, labeller = label_both) +
    #   ylim(0, 0.3) +
    #   scale_x_continuous(breaks = unique(dset_es_cl$test_time)) +
    #   xlab("Age at test") +
    #   ggtitle("Test positive given transition early clinical stage")
    ggplot(data = dset_es_cl, aes(x = m1, y = prosp_sens, color = factor(m2))) +
      # geom_point(aes(size = fraction), alpha = 0.8) +
      geom_line() +
      facet_grid(.~test_time) +
      xlab("MST from early stage to clinical") +
      ylab("Sensitivity in early stage") +
      scale_color_discrete(name = "MST from early stage to late stage") +
      scale_size_continuous(name = "Fraction of early stage to clinical") +
      theme(legend.position = "bottom",
            legend.box = "vertical") +
      scale_y_continuous(limits=c(0, 0.8),
                         breaks=seq(0, 1, by=0.2),
                         labels=label_percent(accuracy=1),
                         expand=c(0, 0))
    if(saveit){
      filename <- str_glue("fig_analytic_multistate_es_cl_{datastamp}.{ext}")
      ggsave(here("figure", filename), width = 6, height = 4)
    }
  }
    
  # Case 2: test positive among early stage to clinical or early stage to late stage
  if(trans_type == 2){
    prosp_sens_es_cl_ls <- foreach(i = cset$i, .combine = "rbind") %dopar% {
      test_es_cl <- test_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i])
      test_es_ls <- test_es_ls_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i])
      case_es_cl <- case_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i])
      case_es_ls <- case_es_ls_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i])
      sens <- (test_es_cl + test_es_ls) / (case_es_cl + case_es_ls)
      sens_es <- test_es_cl/case_es_cl
      sens_ls <- test_es_ls/case_es_ls
      frac <- case_es_cl / (case_es_cl + case_es_ls)
      data.frame(prosp_sens_es = c(sens_es),
                 prosp_sens_ls = c(sens_ls),
                 prosp_sens = c(sens),
                 fraction = c(frac))
    }
    dset_es_cl_ls <- data.frame(cset, prosp_sens_es_cl_ls)
    # initial age is 40
    dset_es_cl_ls <- dset_es_cl_ls %>% mutate(test_time = (40 + test_time))
    # dset_es_cl_ls_55 <- dset_es_cl_ls%>%filter(test_time==55)
    # write.csv(dset_es_cl_ls_55, here("tables", "tbl_sens_multistate_55.csv"))
    
    # ggplot(data = dset_es_cl_ls, aes(x = test_time, y = prosp_sens, color = factor(m1))) +
    #   geom_line() +
    #   facet_grid(. ~ m2, labeller = label_both) +
    #   scale_x_continuous(breaks = unique(dset_es_cl$test_time)) +
    #   ylim(0, 0.8) +
    #   xlab("Age at test") +
    #   ggtitle("Test positive given transition to early clinical or late stage")
    
    gg1 <- ggplot(data = dset_es_cl_ls, aes(x = m1, y = prosp_sens, color = factor(m2))) +
      # geom_point(aes(size = fraction), alpha = 0.8) +
      geom_line() +
      facet_grid(.~test_time) +
      xlab("MST from early stage to clinical") +
      ylab("Sensitivity in early stage") +
      ggtitle("Analytic formula") +
      scale_x_continuous(breaks = m1) +
      scale_color_discrete(name = "MST from early stage to late stage") +
      theme(legend.position = "bottom",
            legend.box = "vertical") +
      scale_y_continuous(limits=c(0, 0.8),
                         breaks=seq(0, 1, by=0.2),
                         labels=label_percent(accuracy=1),
                         expand=c(0, 0))
    
    gg2 <- dset_es_cl_ls %>%
      filter(test_time == 55) %>%
      ggplot(aes(x = m1, y = m2, fill = prosp_sens)) + 
      geom_tile(aes(height = fraction*5, width = fraction*5)) +
      geom_text(aes(label = label_percent(accuracy = 1)(prosp_sens)), color = "white") +
      # geom_text(colour = "white", nudge_y = -1, aes(label = paste0('Proportion: ', round((fraction * 100)), '%'))) +
      scale_fill_gradient(low = "grey70", high = "grey20") +
      scale_x_continuous(breaks = m1, expand = c(0, 0)) + 
      scale_y_continuous(breaks = m2, expand = c(0, 0)) + 
      labs(
        title = "Sensitivity in early preclinical stage \nTest at age 55", 
        x = expression(paste("Mean years from early stage to early clinical (", m[1], ")")), 
        y = expression(paste("Mean years from early to late stage (", m[2], ")"))) +
      theme(legend.position = "none") 


    if(saveit){
      # filename1 <- str_glue("fig_analytic_multistate_es_cl_ls_{datastamp}.{ext}")
      # ggsave(here("figure", filename1), width = 6, height = 4)
      filename2 <- str_glue("fig_analytic_multistate_heatmap_{datastamp}.{ext}")
      ggsave(filename = here("figure", filename2), 
             plot = print(gg2),
             width = 8, height = 8)
    }
  }
  
  # Case 3: test positive among all preclinical
  if(trans_type == 3){
    prosp_sens_all_preclinical <- foreach(i = cset$i, .combine = "rbind") %dopar% {
      (test_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i]) + 
         test_es_ls_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i]) +
         test_ls_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i], m3 = m3)) / 
        (case_es_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i]) + 
           case_es_ls_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i]) +
           case_ls_cl_integral(t = cset$test_time[i], preonset_rate = rate, m1 = cset$m1[i], m2 = cset$m2[i], m3 = m3))
    }
    dset_all_preclinical <- data.frame(cset, prosp_sens = c(prosp_sens_all_preclinical))
    # initial age is 40
    dset_all_preclinical <- dset_all_preclinical %>% mutate(test_time = (40 + test_time))
    
    ggplot(data = dset_all_preclinical, aes(x = m1, y = prosp_sens, color = factor(m2))) +
      geom_line() +
      facet_grid(.~test_time) +
      xlab("MST from early stage to clinical") +
      ylab("Sensitivity in early stage") +
      ggtitle("Analytic formula") +
      scale_x_continuous(breaks = m1) +
      scale_color_discrete(name = "MST from early stage to late stage") +
      theme(legend.position = "bottom",
            legend.box = "vertical") +
      scale_y_continuous(limits=c(0, 0.8),
                         breaks=seq(0, 1, by=0.2),
                         labels=label_percent(accuracy=1),
                         expand=c(0, 0))
    if(saveit){
      filename <- str_glue("fig_analytic_multistate_all_preclinical_{datastamp}.{ext}")
      ggsave(here("figure", filename), width = 6, height = 4)
    }
  }
  
}

control <- function(test_time = seq(10, 20, 5),
                    pre_onset_rate = 0.1,
                    mst_es_cl = c(2, 6, 10),
                    mst_es_ls = c(2, 6, 10),
                    mst_ls_cl = 2,
                    ext = "png",
                    saveit){
  fig_prosp_sens(test_time = test_time, 
                 rate = pre_onset_rate, 
                 m1 = mst_es_cl, 
                 m2 = mst_es_ls, 
                 m3 = mst_ls_cl,
                 trans_type = 2,
                 ext = ext,
                 saveit = saveit)
}
control(saveit = TRUE)
 



