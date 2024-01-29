library(here)
library(foreach)
library(doParallel)

source(here("src", "analytic_sens_multistate_funs.R"))


test_time <- seq(45, 55, 0.5)
rate = 0.2
start = 1
end = 2
lower = 0 
upper = 0.3 
m1 = 4
m2 = 2
m3 = 1
# plot prospective sensitivity for test @ between early stage and clinical onset
fig_prosp_sens <- function(test_time, rate, start, end, lower, upper, m1, m2, m3 = NULL){
  ############
  # test_time: test time at t
  # rate: preclinical onset rate
  # start: initial state
  # end: end state
  # lower: lower bound of true sensitivity
  # upper: upper bound of true sensitivity
  # m1: mean sojourn time for transition from early stage to clinical onset
  # m2: mean sojourn time for transition from early stage to late stage
  # m3: mean sojourn time for transition from late stage to clinical onset
  ############
  no_cores <- detectCores() - 1  # leaving one core free
  registerDoParallel(cores=no_cores)
  
  # early state to clinical onset
  prosp_sens_es_cl <- foreach(t = test_time, .combine = "rbind") %dopar% {
    test_es_cl_integral(t, preonset_rate = rate, m1, m2) / case_es_cl_integral(t, preonset_rate = rate, m1, m2)
  }
  # early state to late stage
  prosp_sens_es_cl_ls <- foreach(t = test_time, .combine = "rbind") %dopar% {
    (test_es_cl_integral(t, preonset_rate = rate, m1, m2) + test_es_ls_integral(t, preonset_rate = rate, m1, m2)) / 
      (case_es_cl_integral(t, preonset_rate = rate, m1, m2) + case_es_ls_integral(t, preonset_rate = rate, m1, m2))
  }
  
}

