## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(learnr)
knitr::opts_chunk$set(echo = FALSE)

# load packages
library(learnr)
library(gradethis)
library(sortable)
library(tidyverse)
library(learnrhash) #devtools::install_github("rundel/learnrhash")
library(showtext)
library(googlesheets4)
library(mcreplicate)
library(knitr)
library(hrbrthemes)
library(here)
library(lme4)
library(progressr)
library(janitor)
library(future)
library(fixest)
library(broom)
lag <- dplyr::lag
#devtools::install_github("graveja0/HPOL8539PKG")
#library(HPOL8539PKG)
# devtools::load_all("../../HPOL8539PKG")


# don't echo chunks
knitr::opts_chunk$set(echo = FALSE)

# apply theme to ggplot
ggplot2::theme_set(theme_bw())

map_multicore <- function(.x, .f, ..., .id = NULL) {
  .f <- purrr::as_mapper(.f, ...)
  p <- progressor(steps = length(.x))
  f <- function(...) {
    p()
    .f(...)
  }
  furrr::future_map(.x, f, ..., .id = .id)
}

plot_sampling_distribution <- function(x,truth) {
  d <- density(x)
  p_df <- as_tibble(cbind(x = d$x, density = d$y))
  p_df %>%
    ggplot(aes(x = x, y = density)) + geom_line() +
    #hrbrthemes::theme_ipsum() +
    labs(x = "Estimate", y = "Density") +
    geom_vline(aes(xintercept = truth)) +
    annotate("text",x = mean(x), y = min(d$y*1.2), vjust=-1,label  = glue::glue("  \tMean: {formatC(mean(x),digits = 3, format='f')}\n   SD: {formatC(sd(x),digits = 3, format = 'f')}"), hjust = 0)
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
    #theme_ipsum() +
    scale_y_discrete(breaks = NULL) +
    geom_vline(aes(xintercept = truth)) +
    labs(title= glue("Confidence Intervals for {prettyNum(length(K),big.mark=',')} of {prettyNum(length(res$m),big.mark=',')} Estimates"),
         y= "Sampling Iteration",x = "Estimate",
         subtitle= glue("{formatC(100*mean(res$covered),digits = 1, format='f')}% of confidence intervals cover the truth")) +
    scale_colour_manual(values = c("black","red")) +
    theme(legend.position = "none")
}

options("scipen" = 100, "digits" = 5)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
params =list(
  N = 1000,
  mean_x_i = 2,
  sd_x_i = 0.5,
  beta = 1,
  tau = 0.5,
  sigma_sq_epsilon = 1
)

params %>% data.frame() %>% gather(param,value) %>% kable(caption = "Parameter Values (Eq. 1)")


## ----basic-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ##################################
## # Step 1: Parameterize the Problem
## ##################################
## params =list(
##   N = 1000,
##   mean_x_i = 2,
##   sd_x_i = 0.5,
##   beta = 1,
##   tau = 0.5,
##   sigma_sq_epsilon = 1
## )
## 
## ####################################################
## # Step 2: Define a Data Generation Process Function
## #####################################################
## dgp_df_u =function(params) {
##   with(params,{
##     df =
##       tibble(
##         # Covariate
##         x_i = rnorm(N, mean_x_i, sd_x_i),
##         # ADD UNOBSERVED HETEROGENEITY TERM
##         u_i = rnorm(N, mean = 0, sd = 1)) %>%
##         # Induce correlation between u_i and treatment;
##         # higher values of u_i make it more likely you're treated.
##         rowwise() %>% # This allows us to get each value's pr_treated in the line below.
##         mutate(pr_treated = boot::inv.logit(u_i)) %>%
##         ungroup() %>%  # This undoes the rowwise
##         # Treatment indicator
##         mutate(d_i = rbinom(N, size = 1, prob = pr_treated)) %>%
##         mutate(epsilon_i = rnorm(N, mean = 0, sd = sigma_sq_epsilon)) %>%
##         # u_i is also in the DGP for y_i
##         #mutate(y_i = beta * x_i + tau * d_i + u_i + epsilon_i) %>%
##         mutate(y_i = beta * x_i + tau * d_i  + epsilon_i) %>%
##         # because u_i is unobserved, we strip it from the "observed" data output.
##         select(-u_i)
##     return(df)
##   })
## }
## 
## ############################################
## # 3. Define and Apply an Estimation Function
## #############################################
## estimator_fn = function(df) {
##   out =
##     df %>%
##       lm(y_i ~ x_i + d_i, data = .)
##   return(out)
## }
## 
## ########################################
## # 4. Define the discriminator function
## # (For this exercise we want to extract tau-hat
## # i.e., the coefficient on treated)
## ########################################
## disc_fn = function(fit) {
##   fit_ =broom::tidy(fit)   # This cleans up the fitted regression object
##   out =fit_ %>%
##     filter(term=="d_i") %>%
##     pull(estimate)
## 
##   return(out)
## }
## 
## ###############################################
## # 5. Define a compound function that executes
## # steps 1-4 based on the parameter inputs.
## ###############################################
## 
## generate_estimate_discriminate <- function(params) {
##   params %>% # Step 1: Parameterize the problem
##       dgp_df_u() %>%  # Step 2: Define the data generation process
##         estimator_fn() %>%  # Step 3: Estimate
##           disc_fn() %>% # Step 4: Pull out what you need
##             data.frame(tau_hat = .) # store the result as a data frame object
## }
## 
## # Monte Carlo simulation based on 100 different realizations of the DGP:
## M = 100
## set.seed(123)
## mc_result <- 1:M %>% map_df(~generate_estimate_discriminate(params))
## 
## plot_sampling_distribution(mc_result$tau_hat, truth = params$tau)
## ggsave(here::here("Panel Data Methods I/images/samp-dist-noUi.png"),width=6,height=5)
## 
## 


## ----out.width="50%", echo = FALSE, fig.align="center",fig.cap="Sampling Distribution of OLS Estimator of tau"---------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/samp-dist-noUi.png")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## params %>%
##   generate_data() %>%
##     estimate() %>%
##       discriminate()
## 
## do_it_all <- function(params) {
##   params %>%
##     generate_data() %>%
##       estimate() %>%
##         discriminate()
## }
## 


## ----ovb1--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##################################
# Step 1: Parameterize the Problem
##################################
params =list(
  N = 1000,
  mean_x_i = 2,
  sd_x_i = 0.5,
  beta = 1,
  tau = 0.5,
  sigma_sq_epsilon = 1
)

####################################################
# Step 2: Define a Data Generation Process Function
#####################################################
dgp_df_u =function(params) {
  with(params,{
    df =
      tibble(
        # Covariate
        x_i = rnorm(N, mean_x_i, sd_x_i),
        # ADD UNOBSERVED HETEROGENEITY TERM
        u_i = rnorm(N, mean = 0, sd = 1)) %>% 
        # Induce correlation between u_i and treatment;
        # higher values of u_i make it more likely you're treated. 
        rowwise() %>% # This allows us to get each value's pr_treated in the line below. 
        mutate(pr_treated = boot::inv.logit(u_i)) %>% 
        ungroup() %>%  # This undoes the rowwise 
        # Treatment indicator
        mutate(d_i = rbinom(N, size = 1, prob = pr_treated)) %>% 
        mutate(epsilon_i = rnorm(N, mean = 0, sd = sigma_sq_epsilon)) %>% 
        # u_i is also in the DGP for y_i 
        mutate(y_i = beta * x_i + tau * d_i + u_i + epsilon_i) %>% 
        # because u_i is unobserved, we strip it from the "observed" data output. 
        select(-u_i)
    return(df)
  })
}

############################################
# 3. Define and Apply an Estimation Function
#############################################
estimator_fn = function(df) {
  out =
    df %>% 
      lm(y_i ~ x_i + d_i, data = .)
  return(out)
}

########################################
# 4. Define the discriminator function 
# (For this exercise we want to extract tau-hat 
# i.e., the coefficient on treated)
########################################
disc_fn = function(fit) {
  fit_ =broom::tidy(fit)   # This cleans up the fitted regression object
  out =fit_ %>% 
    filter(term=="d_i") %>% 
    pull(estimate)
  
  return(out)
}

###############################################
# 5. Define a compound function that executes
# steps 1-4 based on the parameter inputs. 
###############################################

generate_estimate_discriminate <- function(params) {
  params %>% # Step 1: Parameterize the problem
      dgp_df_u() %>%  # Step 2: Define the data generation process
        estimator_fn() %>%  # Step 3: Estimate 
          disc_fn() %>% # Step 4: Pull out what you need
            data.frame(tau_hat = .) # store the result as a data frame object
}




## ----mc_ovb, exercise.setup = "ovb1"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Monte Carlo simulation based on 100 different realizations of the DGP:
M = 100
result_lm <- 1:M %>% map_df(~generate_estimate_discriminate(params))

plot_sampling_distribution(result_lm$tau_hat, truth = params$tau)


## ----dgp_panel_setup---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
params_panel <- list(
  N = 1000,
  T = 2,
  tx_time = 2, 
  rho_t = 0.8,
  beta_0 = 0.5,
  beta_1 = 2,
  tau = 0.5,
  p_d = 0.5
)

dgp_panel <- function(params) {
  with(params, {

    # Time effects
    t_ <-
      data.frame(t = 1:T,
                 gamma_t = arima.sim(n=T, list(ar = rho_t, order=c(1,0,0))) %>% as.vector())

    # Individual measures and effects
    i_ <-
      data.frame(
        unit_id = 1:N,
        x_i = rnorm(N, mean = 0, sd = 1),
        u_i = rnorm(N, mean = 0, sd = 1)) %>%
      rowwise() %>% # This allows us to get each value's pr_treated in the line below. 
      mutate(pr_treated = boot::inv.logit(u_i)) %>% 
      ungroup() %>%  # This undoes the rowwise 
      # Treatment indicator
      mutate(d_i = rbinom(N, size = 1, prob = pr_treated)) %>% 
      ungroup()

    crossing(unit_id = i_$unit_id,t = t_$t) %>%
      left_join(i_,"unit_id") %>%
      left_join(t_,"t") %>%
      mutate(d_i = ifelse(t<tx_time,0,d_i)) %>%
      mutate(y_i = beta_0 + beta_1 * x_i + tau * d_i + u_i + gamma_t + rnorm(N, mean = 0, sd = 1))
  })
}

params_panel %>% 
  dgp_panel()



## ----pooledOLS1, exercise.setup="dgp_panel_setup"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimator_fn_pols <- function(df) {
  lm(y_i ~ x_i + d_i, data = df)
}

params_panel %>% 
  dgp_panel() %>% 
    estimator_fn_pols()



## ----pooledOLS2, exercise.setup="pooledOLS1"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define a discrimination function for the POLS estimator 
disc_fn_pols = function(fit) {
  fit_ =broom::tidy(fit)   
  out =fit_ %>% 
    filter(term=="d_i") %>% 
    pull(estimate)
  return(out)
}

generate_estimate_discriminate_pols <- function(params) {
  params %>% # Step 1: Parameterize the problem
      dgp_panel() %>%  # Step 2: DGP
        estimator_fn_pols() %>%  # Step 3: Estimate 
          disc_fn_pols() %>% # Step 4: Coefficient of interest
            data.frame(tau_hat = .) # Store as data fram
}

M = 100
result_pols <- 1:M %>% map_df(~generate_estimate_discriminate_pols(params_panel))
plot_sampling_distribution(result_pols$tau_hat, truth = params_panel$tau)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## dgp_panel_randomizedD <- function(params) {
##   with(params, {
## 
##     # Time effects
##     t_ <-
##       data.frame(t = 1:T,
##                  gamma_t = rnorm(n=T))
##                  #gamma_t = arima.sim(n=T, list(ar = rho_t, order=c(1,0,0))) %>% as.vector())
## 
##     # Individual measures and effects
##     i_ <-
##       data.frame(
##         unit_id = 1:N,
##         x_i = rnorm(N, mean = 0, sd = 1),
##         u_i = rnorm(N, mean = 0, sd = 1)) %>%
##       rowwise() %>% # This allows us to get each value's pr_treated in the line below.
##       mutate(pr_treated = 0.5) %>%
##       ungroup() %>%  # This undoes the rowwise
##       # Treatment indicator
##       mutate(d_i = rbinom(N, size = 1, prob = pr_treated)) %>%
##       ungroup()
## 
##     crossing(unit_id = i_$unit_id,t = t_$t) %>%
##       left_join(i_,"unit_id") %>%
##       left_join(t_,"t") %>%
##       mutate(d_i = ifelse(t<tx_time,0,d_i)) %>%
##       mutate(y_i = beta_0 + beta_1 * x_i + tau * d_i + u_i + gamma_t + rnorm(N, mean = 0, sd = 1))
##   })
## }
## 
## generate_estimate_discriminate_pols_randomizedD <- function(params) {
##   params %>% # Step 1: Parameterize the problem
##       dgp_panel_randomizedD() %>%  # Step 2: Define the data generation process
##         estimator_fn_pols() %>%  # Step 3: Estimate
##           disc_fn_pols() %>% # Step 4: Pull out what you need
##             data.frame(tau_hat = .) # store the result as a data frame object
## }
## 
## M = 100
## result_pols_randomizedD <- 1:M %>% map_df(~generate_estimate_discriminate_pols_randomizedD(modifyList(params_panel,list(N=10000))))
## plot_sampling_distribution(result_pols_randomizedD$tau_hat, truth = params_panel$tau)
## 


## ----re1, exercise.setup="dgp_panel_setup"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

estimator_fn_re <- function(df) {
   lmer(y_i ~ x_i+ d_i + (1|unit_id) + (1|t), df)
}

params_panel %>% 
  dgp_panel() %>% 
    estimator_fn_re()



## ----re2, exercise.setup="re1"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define a discriminator function
disc_fn_re <- function(fit) {
  fit %>% summary() %>% pluck("coefficients") %>%
    data.frame() %>%
    rownames_to_column() %>%
    janitor::clean_names() %>%
    filter(rowname=="d_i") %>%
    pull(estimate) %>%
    as.vector()
}

# Bundle it all together in one function. 
generate_estimate_discriminate_re <- function(params) {
  suppressWarnings({
    suppressMessages({
      params %>% # Step 1: Parameterize the problem
        dgp_panel() %>%  # Step 2: Define the data generation process
          estimator_fn_re() %>%  # Step 3: Estimate 
            disc_fn_re() %>% # Step 4: Pull out what you need
              data.frame(tau_hat = .) # store the result as a data frame object
    })
  })
}

# Run it 100 times!
M = 100
result_re <- 1:M %>% map_df(~generate_estimate_discriminate_re(params_panel))

plot_sampling_distribution(result_re$tau_hat, truth = params_panel$tau)


## ----dummy1, exercise.setup="dgp_panel_setup"--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimator_fn_dummy <- function(df) {
  lm(y_i ~ d_i + factor(t) + factor(unit_id), df)
}

fit_dummy <- 
  params_panel %>% 
    dgp_panel() %>% 
      estimator_fn_dummy()

fit_dummy %>% 
  broom::tidy() 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## generate_estimate_discriminate_dummy <- function(params) {
##   params %>% # Step 1: Parameterize the problem
##       dgp_panel() %>%  # Step 2: Define the data generation process
##         estimator_fn_dummy() %>%  # Step 3: Estimate
##           disc_fn() %>% # Step 4: Pull out what you need
##             data.frame(tau_hat = .) # store the result as a data frame object
## }
## 
## M = 1000
## result_dummy <-
##   1:M %>%
##   map_multicore(~generate_estimate_discriminate_dummy(params_panel)) %>%
##   bind_rows()
## 
## plot_sampling_distribution(result_dummy$tau_hat, truth = params$tau)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## #ggsave(here::here("Panel Data Methods I/images/samp-dist-dummy.png"),width=6,height=5)


## ----out.width="50%", echo = FALSE, fig.align="center"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/samp-dist-dummy.png")


## ----fistdiff1, exercise.setup="dgp_panel_setup"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
construct_first_diff <- function(df) {
  df_ <- 
    df %>% 
      arrange(unit_id,t) %>% 
      group_by(unit_id) %>% 
      mutate(y_fd = y_i - lag(y_i),
             x_fd = x_i - lag(x_i),
             d_fd = d_i - lag(d_i)) %>% 
      filter(t==2)
  return(df_)
}

estimate_first_diff <- function(df) {
  lm(y_fd ~ d_fd , data = df)
}

params_panel %>% 
  dgp_panel() %>% 
    construct_first_diff() %>% 
      estimate_first_diff() %>% 
        summary() 



## ----fistdiff2, exercise.setup="fistdiff1"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
disc_fn_first_diff = function(fit) {
  fit_ =broom::tidy(fit)   # This cleans up the fitted regression object
  out =fit_ %>% 
    filter(term=="d_fd") %>% 
    pull(estimate)
  
  return(out)
}

generate_estimate_discriminate_first_diff <- function(params) {
  params_panel %>% 
    dgp_panel() %>% 
      construct_first_diff() %>% 
        estimate_first_diff() %>% 
          disc_fn_first_diff() %>% 
              data.frame(tau_hat = .) # store the result as a data frame object
}

 
M = 100
result_first_diff <- 
  1:M %>% 
  map_df(~generate_estimate_discriminate_first_diff(params_panel)) 

plot_sampling_distribution(result_first_diff$tau_hat, truth = params_panel$tau)



## ----demean1,exercise.setup="dgp_panel_setup"--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
demean_data <- function(df) {
  df_ <- 
    df %>% 
      mutate(
            # Step 1: Add in the global mean.
            y_dm = y_i + mean(y_i),
            d_dm = d_i + mean(d_i)) %>% 
            # Step 2: Subtract out the time-period means
            group_by(t) %>% 
            mutate(y_dm = y_dm - mean(y_dm),
            d_dm = d_dm - mean(d_i)) %>% 
            # Step 3: Subtract out the unit-level means. 
            group_by(unit_id) %>% 
            mutate(y_dm = y_dm - mean(y_dm),
            d_dm = d_dm - mean(d_dm)) %>% 
            ungroup()

  return(df_)
}

params_panel %>% 
  dgp_panel() %>% 
    demean_data() 


## ----demean2,exercise.setup="demean1"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimator_fn_dm <- function(df) {
  lm(y_dm ~ d_dm   , data = df)
}

set.seed(123)
params_panel %>% 
  dgp_panel() %>% 
    demean_data() %>% 
      estimator_fn_dm()


## ----demean3,exercise.setup="demean2"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimator_fn_dm2 <- function(df) {
  feglm(y_i ~  d_i | t + unit_id, df, family = "gaussian")
}

set.seed(123)
params_panel %>% 
  dgp_panel() %>% 
    demean_data() %>% 
      estimator_fn_dm2()


## ----demean4,exercise.setup="demean3"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define a discriminator function that collects estimates of \hat \tau
disc_fn_dm = function(fit) {
  fit_ =broom::tidy(fit)   # This cleans up the fitted regression object
  out =fit_ %>% 
    filter(term=="d_dm") %>% 
    pull(estimate)
  
  return(out)
}

# Bundle it all together in a single function. 
generate_estimate_discriminate_dm <- function(params) {
  params %>% 
    dgp_panel() %>% 
        demean_data() %>% 
            estimator_fn_dm() %>% 
              disc_fn_dm() %>% 
                data.frame(tau_hat = .) # store the result as a data frame object
}

M = 100
result_dm <- 
  1:M %>% 
  map_df(~generate_estimate_discriminate_dm(params_panel)) 

plot_sampling_distribution(result_dm$tau_hat, truth = params_panel$tau)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# incidental parameters: https://stats.stackexchange.com/questions/185998/incidental-parameter-problem


## ----cre1, exercise.setup="dgp_panel_setup"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prepare_cre <- function(df) {
  df %>% 
    group_by(t) %>%
    mutate(d_bar_t = mean(d_i),
           x_bar_t = mean(x_i)) %>%
    group_by(unit_id) %>%
    mutate(d_bar_i = mean(d_i),
           x_bar_i = mean(x_i))
}

estimator_fn_cre <- function(df) {
  suppressWarnings({suppressMessages({
    lmer(y_i ~ d_i + d_bar_t + d_bar_i + (1|t) + (1|unit_id), df)
  })})
}

set.seed(123)
params_panel %>% 
  dgp_panel()  %>% 
    prepare_cre() %>% 
      estimator_fn_cre()



## ----cre2, exercise.setup="cre1"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimator_fn_fe <- function(df) {
  suppressWarnings({suppressMessages({
    feglm(y_i ~  d_i | t + unit_id, df, family = "gaussian")
  })})
}

set.seed(123)
params_panel %>% 
  dgp_panel()  %>% 
    prepare_cre() %>% 
      estimator_fn_fe()



## ----cre3, exercise.setup="cre2"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
disc_fn_cre <- function(fit) {
  fit %>% summary() %>% pluck("coefficients") %>%
    data.frame() %>%
    rownames_to_column() %>%
    janitor::clean_names() %>%
    filter(rowname=="d_i") %>%
    pull(estimate) %>%
    as.vector()
}

generate_estimate_discriminate_cre <- function(params) {
  params %>% # Step 1: Parameterize the problem
    dgp_panel() %>%  # Step 2: Define the data generation process
      prepare_cre() %>% # Step 2.5: Get the unit- and time-specific means 
      estimator_fn_cre() %>%  # Step 3: Estimate 
        disc_fn_cre() %>% # Step 4: Pull out what you need
        data.frame(tau_hat = .) # store the result as a data frame object
}


M = 100
result_cre <- 1:M %>% map_df(~generate_estimate_discriminate_cre(params_panel))
plot_sampling_distribution(result_cre$tau_hat, truth = params_panel$tau)




