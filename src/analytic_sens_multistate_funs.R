# function of preclinical onset
f <- function(u, preonset_rate){
  preonset_rate*exp(-preonset_rate*u)
}

# function of true sensitivity in early stage
h <- function(t,
              sojourn_time,
              lower,
              upper){
  # t: time from preclincial onset
  beta0 <- lower
  beta1 <- (upper - lower)/sojourn_time
  true_sens <- beta0 + beta1 * t
  true_sens <- pmax(lower, pmin(upper, true_sens))
  return(true_sens)
}


# function of sojourn time
g <- function(s, mean_sojourn_time){
  1/mean_sojourn_time*exp(-1/mean_sojourn_time*s)
}

#########
# Function of mean sojourn time from early stage to clinical stage given not going to late stage
#########
mst_es_cl_f <- function(s1, s2, m1, m2){
  # t0: time at early stage
  # s1: sojourn time from early stage to clinical onset
  # s2: sojourn time from early stage to late stage
  # m1: mst from early stage to clinical stage
  # m2: mst from early stage to late stage
  s1 * g(s = s1, mean_sojourn_time = m1) * g(s = s2, mean_sojourn_time = m2)
}

mst_es_cl_integral_s1 <- function(s2, m1, m2) {
  integrate(f = mst_es_cl_f, 
            lower = 0, upper = s2, 
           s2 = s2,
           m1 = m1, m2 = m2)$value
}
# Compute the nested integral
mst_es_cl_integral_s2 <- function(m1, m2){
  integrate(f = Vectorize(mst_es_cl_integral_s1), 
                 lower = 0, upper = Inf,
                 m1 = m1, m2 = m2)$value
  }

# mst_es_cl_integral_t2(t0, m1, m2)

#########
# Function of mean sojourn time from early stage to late stage given not going to clinical onset
#########
mst_es_ls_f <- function(s1, s2, m1, m2){
  # t0: time at early stage
  # s1: sojourn time from early stage to clinical onset
  # s2: sojourn time from early stage to late stage
  # m1: mst from early stage to clinical stage
  # m2: mst from early stage to late stage
  s2 * g(s = s1, mean_sojourn_time = m1) * g(s = s2, mean_sojourn_time = m2)
}

mst_es_ls_integral_s2 <- function(s1, m1, m2) {
  integrate(f = mst_es_ls_f,
            lower = 0, upper = s1, 
            s1 = s1,
            m1 = m1, m2 = m2)$value
}
# Compute the nested integral
mst_es_ls_integral_s1 <- function(m1, m2){
  integrate(f = Vectorize(mst_es_ls_integral_s2), 
            lower = 0, upper = Inf,
            m1 = m1, m2 = m2)$value
}

# mst_es_ls_integral_t1(t0, m1, m2)

###############
# Analytic function of stage-specific prospective sensitivity at time t
###############
# test between preclinical and preclinical onset
f_es_cl <- function(t, t0, s1, s2, m1, m2, test = TRUE){
  # t0: time at early stage
  # s1: sojourn time from early stage to clinical onset
  # s2: sojourn time from early stage to late stage
  # m1: mst from early stage to clinical stage
  # m2: mst from early stage to late stage
  f_joint <- g(s = s1, mean_sojourn_time = mst_es_cl_integral_s2(m1, m2)) *
    g(s = s2, mean_sojourn_time = mst_es_ls_integral_s1(m1, m2))
  # f_joint <- g(s = s1, mean_sojourn_time = m1) * 
  #   g(s = s2, mean_sojourn_time = m2)
  if(test){
    f_joint <- f_joint * h(t = (t - t0), sojourn_time = s1, lower = 0, upper = 0.3)
  }
  return(f_joint)
}
# integral over s1
test_es_cl_integral_s1 <- function(t, t0, s2, m1, m2){
  integrate(f = f_es_cl, 
            lower = t-t0, upper = s2, 
            t = t, t0 = t0, s2 = s2,
            m1 = m1, m2 = m2)$value
}
# integral over s2
test_es_cl_integral_s2 <- function(t, t0, m1, m2){
  integrate(f = Vectorize(test_es_cl_integral_s1), 
            lower = t-t0, upper = Inf,
            t = t, t0 = t0,
            m1 = m1, m2 = m2)$value
}
# integral over t0
f_es_cl_integral_t0 <- function(t, t0, preonset_rate, m1, m2){
  f(u = t0, preonset_rate) * test_es_cl_integral_s2(t, t0, m1, m2)
}

test_es_cl_integral <- function(t, preonset_rate, m1, m2){
  integrate(f = Vectorize(f_es_cl_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2)$value
}
# test_es_cl_integral(t, preonset_rate, m1, m2)


# test between preclinical and late stage
f_es_ls <- function(t, t0, t1, t2, m1, m2, test = TRUE){
  # t0: time at early stage
  # t1: time at early stage clinical onset
  # t2: time at late stage
  # m1: mst from early stage to clinical stage
  # m2: mst from early stage to late stage
  f_joint <- g(s = (t1 - t0), mean_sojourn_time = mst_es_cl_integral_t2(t0, m1, m2)) * 
    g(s = (t2 - t0), mean_sojourn_time = mst_es_ls_integral_t1(t0, m1, m2))
  if(test){
    f_joint <- f_joint * h(t = (t - t0), sojourn_time = (t2 - t0), lower = 0, upper = 0.8)
  }
  return(f_joint)
}
# integral over t2
test_es_ls_integral_t2 <- function(t, t0, t1, m1, m2){
  integrate(f = f_es_ls, 
            lower = t, upper = t1, 
            t = t, t0 = t0, t1 = t1,
            m1 = m1, m2 = m2)$value
}
# integral over t1
test_es_ls_integral_t1 <- function(t, t0, m1, m2){
  integrate(f = Vectorize(test_es_ls_integral_t2), 
            lower = t, upper = Inf,
            t = t, t0 = t0,
            m1 = m1, m2 = m2)$value
}
# integral over t0
f_es_ls_integral_t0 <- function(t, t0, preonset_rate, m1, m2){
  f(u = t0, preonset_rate) * test_es_ls_integral_t1(t, t0, m1, m2)
}

test_es_ls_integral <- function(t, preonset_rate, m1, m2){
  integrate(f = Vectorize(f_es_ls_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2)$value
}
# test_es_ls_integral(t, preonset_rate, m1, m2)


# test between late stage and clinical stage
f_ls_cl <- function(t, t0, t1, t2, t3, m1, m2, m3, test = TRUE){
  # t0: time at early stage
  # t1: time at early stage clinical onset
  # t2: time at late stage
  # m1: mst from early stage to clinical stage
  # m2: mst from early stage to late stage
  f_joint <- g(s = (t1 - t0), mean_sojourn_time = mst_es_cl_integral_t2(t0, m1, m2)) * 
    g(s = (t2 - t0), mean_sojourn_time = mst_es_ls_integral_t1(t0, m1, m2)) *
    g(s = (t3 - t2), mean_sojourn_time = m3)
  if(test){
    f_joint <- f_joint * 0.8
  }
  return(f_joint)
  
}
# integral over t1
test_ls_cl_integral_t1 <- function(t, t0, t2, t3, m1, m2, m3){
  integrate(f = f_ls_cl, 
            lower = t2, upper = Inf, 
            t = t, t0 = t0, t2 = t2, t3 = t3,
            m1 = m1, m2 = m2, m3 = m3)$value
}
# integral over t2
test_ls_cl_integral_t2 <- function(t, t0, t3, m1, m2, m3){
  integrate(f = Vectorize(test_ls_cl_integral_t1), 
            lower = t0, upper = t3,
            t = t, t0 = t0, t3 = t3, 
            m1 = m1, m2 = m2, m3 = m3)$value
}
# integral over t3
test_ls_cl_integral_t3 <- function(t, t0, m1, m2, m3){
  integrate(f = Vectorize(test_ls_cl_integral_t2), 
            lower = t, upper = Inf,
            t = t, t0 = t0, 
            m1 = m1, m2 = m2, m3 = m3)$value
}
# integral over t0
f_ls_cl_integral_t0 <- function(t, t0, preonset_rate, m1, m2, m3){
  f(u = t0, preonset_rate) * test_ls_cl_integral_t3(t, t0, m1, m2, m3)
}

test_ls_cl_integral <- function(t, preonset_rate, m1, m2, m3){
  integrate(f = Vectorize(f_ls_cl_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2, m3 = m3)$value
}
# test_ls_cl_integral(t, preonset_rate, m1, m2, m3)

# cases between preclinical onset and clinical onset
# integral over s1
case_es_cl_integral_s1 <- function(t, t0, s2, m1, m2){
  integrate(f = f_es_cl, 
            lower = t-t0, upper = s2, 
            t = t, t0 = t0, s2 = s2,
            m1 = m1, m2 = m2,
            test = FALSE)$value
}
# integral over s2
case_es_cl_integral_s2 <- function(t, t0, m1, m2){
  integrate(f = Vectorize(case_es_cl_integral_s1), 
            lower = t-t0, upper = Inf,
            t = t, t0 = t0,
            m1 = m1, m2 = m2)$value
}
# integral over t0
f_case_es_cl_integral_t0 <- function(t, t0, preonset_rate, m1, m2){
  f(u = t0, preonset_rate) * case_es_cl_integral_s2(t, t0, m1, m2)
}

case_es_cl_integral <- function(t, preonset_rate, m1, m2){
  integrate(f = Vectorize(f_case_es_cl_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2)$value
}
# case_es_cl_integral(t, preonset_rate, m1, m2)

# case between preclinical and late stage
# integral over t2
case_es_ls_integral_t2 <- function(t, t0, t1, m1, m2){
  integrate(f = f_es_ls, 
            lower = t, upper = t1, 
            t = t, t0 = t0, t1 = t1,
            m1 = m1, m2 = m2,
            test = FALSE)$value
}
# integral over t1
case_es_ls_integral_t1 <- function(t, t0, m1, m2){
  integrate(f = Vectorize(case_es_ls_integral_t2), 
            lower = t, upper = Inf,
            t = t, t0 = t0,
            m1 = m1, m2 = m2)$value
}
# integral over t0
f_case_es_ls_integral_t0 <- function(t, t0, preonset_rate, m1, m2){
  f(u = t0, preonset_rate) * case_es_ls_integral_t1(t, t0, m1, m2)
}

case_es_ls_integral <- function(t, preonset_rate, m1, m2){
  integrate(f = Vectorize(f_case_es_ls_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2)$value
}
# case_es_ls_integral(t, preonset_rate, m1, m2)

# case between late stage and clinical stage
# integral over t1
case_ls_cl_integral_t1 <- function(t, t0, t2, t3, m1, m2, m3){
  integrate(f = f_ls_cl, 
            lower = t2, upper = Inf, 
            t = t, t0 = t0, t2 = t2, t3 = t3,
            m1 = m1, m2 = m2, m3 = m3,
            test = FALSE)$value
}
# integral over t2
case_ls_cl_integral_t2 <- function(t, t0, t3, m1, m2, m3){
  integrate(f = Vectorize(case_ls_cl_integral_t1), 
            lower = t0, upper = t3,
            t = t, t0 = t0, t3 = t3, 
            m1 = m1, m2 = m2, m3 = m3)$value
}
# integral over t3
case_ls_cl_integral_t3 <- function(t, t0, m1, m2, m3){
  integrate(f = Vectorize(case_ls_cl_integral_t2), 
            lower = t, upper = Inf,
            t = t, t0 = t0, 
            m1 = m1, m2 = m2, m3 = m3)$value
}
# integral over t0
f_case_ls_cl_integral_t0 <- function(t, t0, preonset_rate, m1, m2, m3){
  f(u = t0, preonset_rate) * case_ls_cl_integral_t3(t, t0, m1, m2, m3)
}

case_ls_cl_integral <- function(t, preonset_rate, m1, m2, m3){
  integrate(f = Vectorize(f_case_ls_cl_integral_t0), 
            lower = 0, upper = t,
            t = t,
            preonset_rate = preonset_rate,
            m1 = m1, m2 = m2, m3 = m3)$value
}
# case_ls_cl_integral(t, preonset_rate, m1, m2)

