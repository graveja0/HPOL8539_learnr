library(tidyverse)
library(broom)
library(furrr)
library(progressr)
library(tictoc)
library(glue)
library(hrbrthemes)
library(fixest)
library(lme4)
library(lmerTest)
library(janitor)
library(glmnet)

M = 1000

# Extended two-way fixed effects https://arelbundock.com/posts/2021-09-30-extendedtwfe/
# Matrix completion: https://yiqingxu.org/packages/fect/fect.html,
# Matrix completion: https://arelbundock.com/posts/panel_heterogeneous_effects/
# ML DID https://apoorvalal.github.io/posts/21052022_semiparDID.html
# https://diff.healthpolicydatascience.org/#semiparametric


# General Steps
# 1. Generate
# 2. Estimate
# 3. Discriminate

# Course Flow:
# 1. Basic DGP to show sampling distribution, SE, CI coverage, power.
# 2. Omitted variables bias (Tx correlated with unobserved heterogeneity)
# 3. Panel data methods to remove unobserved heterogeneity: fixed, random, correlated random effects.
# 4. Statistical inference issues: Tx assignment at cluster level
# 5. Heterogeneous tx and differential timing
# 6. Nonlinear DID and functional form considerations (i.e., change the outcome )
# 7. Synthetic control methods
# 8. Bounds and principal stratification
# 9. Geographic RD

# Add bias vs. consistency ? https://eranraviv.com/bias-vs-consistency/
# Principal Stratification : https://files.eric.ed.gov/fulltext/EJ1160773.pdf
# https://scholar.harvard.edu/lmiratrix/resources and https://arxiv.org/pdf/1701.03139.pdf

# Geographic Regression Discontinuity: http://www.lukebornn.com/papers/rischard_jasa_2020.pdf and https://github.com/maximerischard/GeoRDD.jl
# and https://www.cambridge.org/core/journals/political-analysis/article/geographic-boundaries-as-regression-discontinuities/2A59F3077F49AD2B908B531F6E458430
# and https://titiunik.mycpanel.princeton.edu/papers/KeeleLorchPassarellaSmallTitiunik2017-AIE.pdf

map_multicore <- function(.x, .f, ..., .id = NULL) {
  .f <- purrr::as_mapper(.f, ...)
  p <- progressor(steps = length(.x))
  f <- function(...) {
    p()
    .f(...)
  }
  furrr::future_map(.x, f, ..., .id = .id)
}
generate_correlated_normals <- function(N,rho) {
  C <-
    matrix(0, nrow = 2, ncol = 2, dimnames = list(c("x_i","u_i"),c("x_i","u_i")))
  diag(C) = 1
  C[2,1] <- C[1,2] <- rho
  L <- chol(C)
  tau <- diag(c(1,1))
  Lambda <- tau %*% t(L)
  Z <- rbind(rnorm(N,mean=0,sd=1),rnorm(N,mean=0,sd=1))
  out <-
    Lambda %*% Z %>% t() %>%
    data.frame() %>%
    as_tibble() %>%
    set_names(colnames(C))
  return(out)
}
plot_sampling_distribution <- function(x,truth) {
  d <- density(x)
  p_df <- as_tibble(cbind(x = d$x, density = d$y))
  p_df %>%
    ggplot(aes(x = x, y = density)) + geom_line() +
    hrbrthemes::theme_ipsum() +
    labs(x = "Estimate", y = "Density") +
    geom_vline(aes(xintercept = truth)) +
    annotate("text",x = mean(x), y = min(d$y*1.1), label  = glue("  \tMean: {formatC(mean(x),digits = 3, format='f')}\n   SD: {formatC(sd(x),digits = 3, format = 'f')}"), hjust = 0)
}
plot_cis <- function(x, K, truth) {
  res <- x %>% bind_rows(.id = "m") %>%
    as_tibble() %>%
    mutate(m = factor(m)) %>%
    mutate(m = fct_reorder(m,estimate, .desc = TRUE)) %>%
    mutate(truth = truth) %>%
    rowwise() %>%
    mutate(covered = as.integer(between(truth,conf.low,conf.high))) %>%
    ungroup() %>%
    mutate(color = ifelse(covered ==1 , "","Rejected"))
  
  K = sample(res$m,100, replace =TRUE)
  res %>%
    filter(m %in% K) %>%
    ggplot() +
    geom_errorbar(aes(xmin =  conf.low, xmax = conf.high, y= m,colour = color)) +
    theme_ipsum() +
    scale_y_discrete(breaks = NULL) +
    geom_vline(aes(xintercept = truth)) +
    labs(title= glue("Confidence Intervals for {prettyNum(length(K),big.mark=',')} of {prettyNum(length(res$m),big.mark=',')} Estimates"),
         y= "Sampling Iteration",x = "Estimate",
         subtitle= glue("{formatC(100*mean(res$covered),digits = 1, format='f')}% of confidence intervals cover the truth")) +
    scale_colour_manual(values = c("black","red")) +
    theme(legend.position = "none")
}

#     plan(multisession,workers = parallel::detectCores()-1)
#     tic()
#     res_parallel <-
#         1:M %>% map_multicore(~{
#             data_generator_function_1(params) %>%
#                 estimator_fn() %>%
#                 discriminator_fn()
#         })
#     toc()

##################################################################
##################################################################
# Session 1: Simulation as a Guide for Study Design and Inference
# DGP: y ~ beta_0 + beta_1 * x_i = beta_2 * w + e
# Notes: x_i uncorrelated with tx
#   x_i ~ N(0,1)
#   tx = binom(p_Tx)
##################################################################
##################################################################

# Step 1: Specify the parameters and define the data generation process
params_1 <- list(
  N = 1000,
  beta_0 = 0.5,
  beta_1 = 2,
  beta_2 = 0,
  p_w = 0.5
)

data_generator_function_1 <- function(params) {
  with(params, {
    data.frame(unit_id = 1:N,
               x_i = rnorm(N, mean = 0, sd = 1),
               w = rbinom(N,size=1, prob = p_w)) %>%
      mutate(y = beta_0 + beta_1 * x_i + beta_2 * w + rnorm(N, mean = 0, sd = 1))
  })
}

# Step 2: Define the estimator and parameter of interest
estimator_fn_1 <- function(df) {
  with(df, {
    lm(y ~ x_i + w)
  })
}

# Do it once!
data_generator_function_1(params_1) %>%
  estimator_fn_1()

# Step 3: What do you want to know?

# a: sampling distribution of the parameter estimate
discriminator_fn_1a <- function(fit) {
  fit %>% broom::tidy() %>% filter(term=="w") %>% pull(estimate)
}

# Do it once!
data_generator_function_1(params_1) %>%
  estimator_fn_1() %>%
  discriminator_fn_1a()

# Do it 1,000 times!

plan(multisession,workers = parallel::detectCores()-1)
with_progress({
  res_1a <-
    1:M %>% map_multicore(~{
      data_generator_function_1(params_1) %>%
        estimator_fn_1() %>%
        discriminator_fn_1a()
    })
})
mean(unlist(res_1a))
plot_sampling_distribution(unlist(res_1a), truth = params_1$beta_2)


# Define a metafunction
generate_estimate_discriminate <- function(params, generator, estimator, discriminator) {
  plan(multisession,workers = parallel::detectCores()-1)
  with_progress({
    res <-
      1:M %>% map_multicore(~{
        generator(params) %>%
          estimator() %>%
          discriminator()
      })
  })
}

params_bias <- list(
  N = 1000,
  mu = 0,
  sigma_sq = .8
)
data_generator_function_bias <- function(params) {
  with(params,{
    data.frame(
      x = rnorm(N,mu,sqrt(sigma_sq))
    ) %>%
      as_tibble()
  })
}

estimator_fn_bias <- function(df) {
  n <- nrow(df)
  
  # df %>%
  #        summarise(var = 1/n * sum((x-mean(x))^2))
  df %>%
    summarise(mu = (1/n) * sum(x) + 10/n)
  
}
discriminator_fn_bias <- function(fit) {
  fit
}

data_generator_function_bias(params_bias) %>%
  estimator_fn_bias() %>%
  discriminator_fn_bias()

sd_N10 <-
  generate_estimate_discriminate(params = modifyList(params_bias,list(N=10)),
                                 generator = data_generator_function_bias,
                                 estimator = estimator_fn_bias,
                                 discriminator = discriminator_fn_bias)
sd_N100 <-
  generate_estimate_discriminate(params = modifyList(params_bias,list(N=100)),
                                 generator = data_generator_function_bias,
                                 estimator = estimator_fn_bias,
                                 discriminator = discriminator_fn_bias)
sd_N1000 <-
  generate_estimate_discriminate(params = modifyList(params_bias,list(N=1000)),
                                 generator = data_generator_function_bias,
                                 estimator = estimator_fn_bias,
                                 discriminator = discriminator_fn_bias)

consistency <-
  c(10,25,50,100,250,500,1000) %>%
  map(~(generate_estimate_discriminate(params = modifyList(params_bias,list(N=.x)),
                                       generator = data_generator_function_bias,
                                       estimator = estimator_fn_bias,
                                       discriminator = discriminator_fn_bias)))
names(consistency) <-  c(10,25,50,100,250,500,1000)
consistency %>%
  map(~(unlist(.x) %>% data.frame() %>% setNames("estimate"))) %>%
  bind_rows(.id = "n") %>%
  mutate(n=as.numeric(n)) %>%
  as_tibble() %>%
  ggplot(aes(x= estimate))  +
  geom_vline(aes(xintercept = 0), lty=3) +
  geom_density() +  labs(x = "Estimate", y = "Density",title = 'N = {closest_state}') + ggthemes::theme_few() +
  transition_states(
    states = n,
    transition_length = 2,
    state_length = 1
  ) +
  enter_fade() +
  exit_shrink() +
  ease_aes('sine-in-out')

# b: Type I Error and Confidence Intervals

discriminator_fn_1b <- function(fit) {
  fit %>% broom::tidy(conf.int =TRUE)  %>%
    filter(term == "w") %>%
    select(estimate, p.value, conf.low, conf.high)
}

res_1b <-
  1:M %>%
  map_multicore(~{
    data_generator_function_1(params_1) %>%
      estimator_fn_1() %>%
      discriminator_fn_1b()
  })

res_1b %>% plot_cis(truth = params_1$beta_2)

# c: Type II Error: Power

discriminator_fn_1c <- function(fit,alpha = 0.05) {
  fit %>% broom::tidy(conf.int =FALSE) %>%
    filter(term=="w") %>%
    
    ### Key information: do we reject at p<0.05?
    ### -> store as a binary
    
    mutate(rejected = as.integer(p.value<alpha))
}

# Change beta_2 truth from 0 to 0.05 for power analysis

params_1c <- modifyList(params_1,
                        list(beta_2 = 0.05))

res_1c <-
  1:M %>%
  map_multicore(~{
    data_generator_function_1(params_1c) %>%
      estimator_fn_1() %>%
      discriminator_fn_1c()
  })

res_1c %>%
  bind_rows() %>%
  summarise(power = mean(rejected))

# We have low power; let's boost the sample size, too.

params_1c_bigN <- modifyList(params_1,
                             list(N=10000,
                                  beta_2 = 0.05))

res_1c_bigN <-
  1:M %>%
  map_multicore(~{
    data_generator_function_1(params_1c_bigN) %>%
      estimator_fn_1() %>%
      discriminator_fn_1c()
  })

res_1c_bigN %>%
  bind_rows() %>%
  summarise(power = mean(rejected))

# Part 2: Omitted Variables Bias

# Step 1: Specify the parameters and define the data generation process
params_1_ovb <- list(
  N = 1000,
  beta_0 = 0.5,
  beta_1 = 2,
  beta_2 = 0,
  p_w = 0.5
)

data_generator_function_1_ovb <- function(params) {
  with(params, {
    data.frame(
      unit_id = 1:N,
      x_i = rnorm(N, mean = 0, sd = 1),
      u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>%
      mutate(pr_w = boot::inv.logit(u_i)) %>%
      mutate(w = rbinom(1,1,prob = pr_w)) %>%
      ungroup() %>%
      mutate(y = beta_0 + beta_1 * x_i + beta_2 * w + u_i + rnorm(N, mean = 0, sd = 1))
  })
}

# Step 2: Define the estimator and parameter of interest
#estimator_fn_1


# Step 3: Define the discriminator function (i.e., the beta_2 estimate)
#discriminator_fn_1a

# Do it once!
data_generator_function_1_ovb(params_1_ovb) %>%
  estimator_fn_1()

# Do it 1,000 times!
M = 1000
res_1_ovb <-
  1:M %>% map_multicore(~{
    data_generator_function_1_ovb(params_1_ovb) %>%
      estimator_fn_1() %>%
      discriminator_fn_1a()
  })
mean(unlist(res_1_ovb))- params_1_ovb$beta_2
plot_sampling_distribution(unlist(res_1_ovb), truth = params_1_ovb$beta_2)


##################################################################
##################################################################
# Session 2: Fixed, Random, and Correlated Random Effects
# DGP: y_it ~ beta_0 + beta_1 * x_i = beta_2 * w_it + u_i + e_it
# Notes: x_i correlated with tx
#   x_i ~ N(0,1)
#   w_it = binom(p_tx), where p_tx = inv.logit(u_i)
##################################################################
##################################################################

# First, we'll define some panel data where the treatment is assigned at the individual level but not correlated
# with unobserved heterogeneity.

params_2 <- list(
  N = 1000,
  T = 2,
  rho_t = 0.8,
  beta_0 = 0.5,
  beta_1 = 2,
  beta_2 = 0,
  p_w = 0.5
)

data_generator_function_2a <- function(params) {
  with(params, {
    
    t_ <-
      data.frame(t = 1:T,
                 gamma_t = arima.sim(n=T, list(ar = rho_t, order=c(1,0,0))) %>% as.vector())
    
    i_ <-
      data.frame(
        unit_id = 1:N,
        x_i = rnorm(N, mean = 0, sd = 1),
        u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>%
      mutate(w = rbinom(1,1,prob = p_w)) %>%
      ungroup()
    
    crossing(unit_id = i_$unit_id,t = t_$t) %>%
      left_join(i_,"unit_id") %>%
      left_join(t_,"t") %>%
      mutate(w = ifelse(t==1,0,w)) %>%
      mutate(y = beta_0 + beta_1 * x_i + beta_2 * w + u_i + gamma_t + rnorm(N, mean = 0, sd = 1))
  })
}

estimator_fn_2fe <- function(df) {
  feglm(y ~  w | t + unit_id, df, family = "gaussian")
}

estimator_fn_2re <- function(df) {
  lmer(y ~ w + (1|t) + (1|unit_id), df)
}

df <- data_generator_function_2a(params_2)
df %>% haven::write_dta("~/Desktop/test.dta")

discriminator_fn_2fe <- function(fit) {
  summary(fit,vcov="twoway") %>% broom::tidy() %>%
    filter(term=="w") %>%
    pull(estimate) %>%
    as.vector()
}

discriminator_fn_2re <- function(fit) {
  fit %>% summary() %>% pluck("coefficients") %>%
    data.frame() %>%
    rownames_to_column() %>%
    janitor::clean_names() %>%
    filter(rowname=="w") %>%
    pull(estimate) %>%
    as.vector()
}

plan(multisession,workers = parallel::detectCores()-1)
res_2afe <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2a(params_2) %>%
      estimator_fn_2fe() %>%
      discriminator_fn_2fe()
  })
(p1 <- res_2afe %>% unlist() %>% plot_sampling_distribution(truth = params_2$beta_2))

res_2are <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2a(params_2) %>%
      estimator_fn_2re() %>%
      discriminator_fn_2re()
  })
res_2are %>% unlist() %>% plot_sampling_distribution(truth = params_2$beta_2) +
  geom_line(
    data = ggplot_build(p1)$data[[1]] %>%
      as_tibble(),
    aes(x = x, y = y, colour = "red")
  ) + theme(legend.position = "none")

# b. Correlated Tx and Unobserved Heterogeneity

data_generator_function_2b <- function(params) {
  with(params, {
    t_ <-
      data.frame(t = 1:T,
                 gamma_t = arima.sim(n=T, list(ar = rho_t, order=c(1,0,0))) %>% as.vector())
    
    
    i_ <-
      data.frame(
        unit_id = 1:N,
        x_i = rnorm(N, mean = 0, sd = 1),
        u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>%
      ################################################
    # mutate(tx = rbinom(1,1,p_Tx)) %>% # this is now correlated with u_i, the unobserved (individual-level) heterogeneity
    mutate(pr_w_fn_u = boot::inv.logit(u_i)) %>%
      mutate(w = rbinom(1,1,prob = pr_w_fn_u)) %>%
      ################################################
    ungroup()
    
    crossing(unit_id = i_$unit_id,t = t_$t) %>%
      left_join(i_,"unit_id") %>%
      left_join(t_,"t") %>%
      mutate(w = ifelse(t==1,0,w)) %>%
      mutate(y = beta_0 + beta_1 * x_i + beta_2 * w + u_i + gamma_t + rnorm(N, mean = 0, sd = 1))
  })
}

plan(multisession,workers = parallel::detectCores()-1)
res_2bfe <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2b(params_2) %>%
      estimator_fn_2fe() %>%
      discriminator_fn_2fe()
  })
(p2 <- res_2bfe %>% unlist() %>% plot_sampling_distribution(truth = params_2$beta_2))

res_2bre <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2b(params_2) %>%
      estimator_fn_2re() %>%
      discriminator_fn_2re()
  })
res_2bre %>% unlist() %>% plot_sampling_distribution(truth = params_2$beta_2) +
  geom_line(
    data = ggplot_build(p2)$data[[1]] %>%
      as_tibble(),
    aes(x = x, y = y, colour = "red")
  ) + theme(legend.position = "none")

estimator_fn_2cre <- function(df) {
  df_ <-
    df %>%
    group_by(t) %>%
    mutate(w_t = mean(w)) %>%
    group_by(unit_id) %>%
    mutate(w_i = mean(w))
  lmer(y ~ w + w_t + w_i + (1|t) + (1|unit_id), df_)
}

res_2bcre <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2b(params_2) %>%
      estimator_fn_2cre() %>%
      discriminator_fn_2re()
  })

res_2bcre %>% unlist() %>% plot_sampling_distribution(truth = params_2$beta_2) +
  geom_line(
    data = ggplot_build(p2)$data[[1]] %>%
      as_tibble(),
    aes(x = x, y = y, colour = "red")
  ) + theme(legend.position = "none")

# Plot 95% CI Coverage

df <- data_generator_function_2b(params_2)

fit <- feols(y ~  w | t + unit_id, df)
summary(fit, vcov="twoway") %>% broom::tidy()
summary(fit, vcov=~unit_id) %>% broom::tidy()
summary(fit, cluster=~unit_id) %>% broom::tidy()
summary(fit) %>% broom::tidy()

discriminator_fn_2fe_type1 <- function(fit, alpha = 0.05) {
  summary(fit,vcov=~unit_id) %>% broom::tidy()  %>%
    filter(term == "w") %>%
    pull(p.value) %>%
    as.vector() -> p
  as.integer(p<alpha)
}

fit <- df %>% estimator_fn_2re()

discriminator_fn_2re_type1 <- function(fit, alpha = 0.05) {
  summary(fit,vcov=~unit_id) %>% broom::tidy()  %>%
    filter(term == "w") %>%
    pull(p.value) %>%
    as.vector() -> p
  as.integer(p<alpha)
}

res_2fe_type1 <-
  1:M %>%
  map_multicore(~{
    data_generator_function_2a(params_2) %>%
      estimator_fn_2fe() %>%
      discriminator_fn_2fe_type1()
  })

mean(unlist(res_2fe_type1))

##################################################################
##################################################################
# Session 3: Statistical Inference for Policy Evaluation
##################################################################
##################################################################

params_3 <- list(
  N = 1000,                # Sample size
  C = 10,                  # Total number of clusters
  N_treated = 5,           # Number of treated clusters
  T = 2008:2012,           # Analytic time window (calendar years)
  
  rho_t = 0.8,             # Global time series correlation parameter
  rho_c = 0.8,             # Cluster-level time series correlation parameter
  
  beta_0 = 0.5,            # Outcome model intercept term
  beta_1 = 2,              # Outcome model coefficient on individual characteristic x_i
  beta_w = list(           # Treatment effect parameters
    "2010" = c(0,0)
  )
)

data_generator_function_3a <- function(params) {
  
  # Treatment is simply a function of unobserved group heterogeneity.
  
  with(params, {
    t_ <-
      data.frame(t = T,
                 gamma_t = arima.sim(n=length(T), list(ar = rho_t, order=c(1,0,0))) %>% as.vector()) %>%
      as_tibble(); t_
    
    treatment_times <-names(beta_w)
    tx_t_prob <- rep(1/length(treatment_times),length(treatment_times))
    tx_t_prob <- tx_t_prob / sum(tx_t_prob)
    c_ <-
      data.frame(c = 1:C,
                 u_c = rnorm(C,mean = 0, sd = 1)) %>%
      mutate(pr_w = boot::inv.logit(u_c)) %>%
      arrange(desc(pr_w)) %>%
      # The top N_treated clusters are considered "ever treated"
      mutate(ever_treated = as.integer(row_number() %in% 1:N_treated)) %>%
      mutate(pop_size =rnbinom(C,size=1,mu=1000)) %>%
      mutate(pr_c = 1/C) %>%
      mutate(c = paste0(c)) %>%
      as_tibble() %>%
      arrange(as.numeric(c)) %>%
      rowwise() %>%
      mutate(treatment_time = ifelse(ever_treated==1,as.numeric(paste0(sample(treatment_times,1,prob=tx_t_prob))),NA)) %>%
      select(c, ever_treated,treatment_time,u_c, pr_w, everything()); c_
    
    pr_c <- c_$pr_c %>% set_names(c_$c)
    
    i_ <-
      data.frame(
        unit_id = 1:N,
        x_i = rnorm(N, mean = 0, sd = 1),
        u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>%
      mutate(c = sample(x= names(pr_c),size = 1, prob = pr_c)) %>%
      left_join(c_,"c") %>%
      ungroup() %>%
      as_tibble(); i_
    
    crossing(unit_id = i_$unit_id,t = t_$t) %>%
      left_join(i_,"unit_id") %>%
      left_join(t_,"t") %>%
      mutate(w = case_when(t >= treatment_time ~ 1, TRUE ~ 0)) %>%
      mutate(time_since_tx = case_when(w==1 ~ t-treatment_time, TRUE~0)) %>%
      mutate(e_it = rnorm(n=nrow(.))) %>%
      rowwise() %>%
      mutate(beta_0_ = beta_0,
             beta_1_ = beta_1,
             beta_2_ = ifelse(ever_treated==1,beta_w[[paste0(treatment_time)]][1],0),
             beta_3_ = ifelse(ever_treated==1,beta_w[[paste0(treatment_time)]][2],0)) %>%
      ungroup() %>%
      mutate(y = beta_0_ + beta_1_ * x_i + beta_2_ * w + beta_3_ * time_since_tx + u_i + u_c + e_it)
    
  })
}

data_generator_function_3a(params_3)

estimator_fn_3fe <- function(df) {
  feglm(y ~  w + factor(t) |  unit_id + c, df, family = "gaussian")
}

data_generator_function_3a(params_3) %>% estimator_fn_3fe()

estimator_fn_3re <- function(df) {
  lmer(y ~ w + (1|t) + (1|unit_id) + (1|c), df)
}

data_generator_function_3a(params_3) %>% estimator_fn_3re()

estimator_fn_3cre <- function(df) {
  df_ <-
    df %>%
    group_by(t) %>%
    mutate(w_t = mean(w)) %>%
    group_by(unit_id) %>%
    mutate(w_i = mean(w))
  lmer(y ~ w + w_t + w_i  + (1|t) + (1|unit_id) + (1|c), df_)
}
data_generator_function_3a(params_3) %>% estimator_fn_3cre()

fit <- data_generator_function_3a(params_3) %>% estimator_fn_3fe()

discriminator_fn_3fe <- function(fit, alpha = 0.05) {
  summary(fit,vcov="twoway") %>% broom::tidy() %>%
    filter(term=="w") %>%
    mutate(reject = as.integer(p.value < alpha)) %>%
    mutate(se = "twoway") %>%
    bind_rows(
      summary(fit,cluster = "c") %>% broom::tidy() %>%
        filter(term=="w") %>%
        mutate(reject = as.integer(p.value < alpha)) %>%
        mutate(se = "cluster")
    )
}

discriminator_fn_3re <- function(fit, alpha =0.05) {
  fit %>% summary() %>% pluck("coefficients") %>%
    data.frame() %>%
    rownames_to_column() %>%
    janitor::clean_names() %>%
    filter(rowname=="w") %>%
    select(term = rowname, estimate,std.error = std_error, statistic = t_value, p.value = pr_t) %>%
    mutate(reject = as.integer(p.value<alpha))
}

plan(multisession,workers = parallel::detectCores()-1)
with_progress({
  res_3afe <-
    1:M %>%
    map_multicore(~{
      data_generator_function_3a(params_3) %>%
        estimator_fn_3fe() %>%
        discriminator_fn_3fe()
    })
  system("say All done!")
})

fdr <- res_3afe %>% bind_rows() %>% filter(se=="cluster") %>% as_tibble() %>% summarise(reject = mean(reject)) %>% pull(reject) %>% {glue("False Discovery Rate = {.*100}%\nalpha = 5%")}

(p3 <- res_3afe %>% bind_rows() %>% filter(se=="cluster") %>% pull(estimate) %>% plot_sampling_distribution(truth = params_3$beta_2) + annotate("text",x=-.25,y = 5, label = fdr,hjust=0,colour = "darkred"))


with_progress({
  res_3are <-
    1:M %>%
    map_multicore(~{
      data_generator_function_3a(params_3) %>%
        estimator_fn_3re() %>%
        discriminator_fn_3re()
    })
  system("say All done!")
})
res_3are %>% bind_rows() %>% pull(estimate) %>% plot_sampling_distribution(truth = params_3$beta_2) +
  geom_line(
    data = ggplot_build(p3)$data[[1]] %>%
      as_tibble(),
    aes(x = x, y = y, colour = "red")
  ) + theme(legend.position = "none")

fdr <- res_3are %>% bind_rows()  %>% as_tibble() %>% summarise(reject = mean(reject)) %>% pull(reject) %>% {glue("False Discovery Rate = {.*100}%\nalpha = 5%")}; fdr

estimator_fn_2cre <- function(df) {
  df_ <-
    df %>%
    group_by(t) %>%
    mutate(w_t = mean(w)) %>%
    group_by(unit_id) %>%
    mutate(w_i = mean(w))
  lmer(y ~ w + w_t + w_i + (1|t) + (1|unit_id) + (1|c), df_)
}

with_progress({
  res_3acre <-
    1:M %>%
    map_multicore(~{
      data_generator_function_3a(params_3) %>%
        estimator_fn_2cre() %>%
        discriminator_fn_3re()
    })
  system("say All done!")
})

fdr <- res_3acre %>% bind_rows()  %>% as_tibble() %>% summarise(reject = mean(reject)) %>% pull(reject) %>% {glue("False Discovery Rate = {.*100}%\nalpha = 5%")}; fdr


res_3acre %>% bind_rows() %>% pull(estimate) %>% plot_sampling_distribution(truth = params_3$beta_2) +
  geom_line(
    data = ggplot_build(p3)$data[[1]] %>%
      as_tibble(),
    aes(x = x, y = y, colour = "red")
  ) + theme(legend.position = "none")


data_generator_function_3b <- function(params) {
  
  # Treatment is simply a function of unobserved group heterogeneity.
  
  with(params, {
    
    n_treated <- p_w * C
    
    t_ <-
      data.frame(t = 1:T,
                 gamma_t = arima.sim(n=T, list(ar = rho_t, order=c(1,0,0))) %>% as.vector())
    
    c_ <-
      data.frame(c = 1:C,
                 u_c = rnorm(C,mean = 0, sd = 1)) %>%
      mutate(pr_x = boot::inv.logit(u_c)) %>%
      arrange(desc(pr_x)) %>%
      mutate(w = as.integer(row_number() %in% 1:n_treated)) %>%
      mutate(pop_size =rnbinom(C,size=1,mu=100)) %>%
      mutate(pr_c = pop_size/sum(pop_size)) %>%
      mutate(c = paste0(c))
    
    ar1_c <-
      c_$c %>%
      map(~(arima.sim(n=T,list(ar = rho_c,order=c(1,0,0))) %>% as.vector() %>% data.frame()))  %>%
      bind_rows(.id = "c") %>%
      group_by(c) %>%
      mutate(t = row_number()) %>%
      set_names(c("c","u_ct","t")) %>%
      mutate(c = paste0(c),
             t = as.numeric(paste0(t)) )
    
    pr_c <- rep(1/C,C) %>% set_names(c_$c)
    
    i_ <-
      data.frame(
        unit_id = 1:N,
        x_i = rnorm(N, mean = 0, sd = 1),
        u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>%
      mutate(c = sample(x= names(pr_c),size = 1, prob = pr_c)) %>%
      left_join(c_,"c") %>%
      ungroup()
    
    crossing(unit_id = i_$unit_id,t = t_$t) %>%
      left_join(i_,"unit_id") %>%
      left_join(t_,"t") %>%
      left_join(ar1_c,c("c","t")) %>%
      mutate(w = ifelse(t==1,0,w)) %>%
      mutate(y = beta_0 + beta_1 * x_i + beta_2 * w + u_i + gamma_t + u_c + u_ct + rnorm(N, mean = 0, sd = 1))
  })
}


plan(multisession,workers = parallel::detectCores()-1)
with_progress({
  res_3bfe <-
    1:M %>%
    map_multicore(~{
      data_generator_function_3b(modifyList(params_3,list(C=4))) %>%
        estimator_fn_3fe() %>%
        discriminator_fn_3fe()
    })
  system("say All done")
})
fdr <- res_3bfe %>% bind_rows() %>% filter(se=="twoway") %>% as_tibble() %>% summarise(reject = mean(reject)) %>% pull(reject) %>% {glue("False Discovery Rate = {.*100}%\nalpha = 5%")}; fdr

(p3 <- res_3bfe %>% bind_rows() %>% filter(se=="cluster") %>% pull(estimate) %>% plot_sampling_distribution(truth = params_3$beta_2) + annotate("text",x=-.25,y = .25, label = fdr,hjust=0,colour = "darkred"))


