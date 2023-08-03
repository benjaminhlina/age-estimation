# ---- load packages -----

library(data.table)
library(dplyr)
library(here)
library(forcats)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(ggmcmc)
library(openxlsx)
library(postpack)
library(patchwork)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# ---- read data -----

lt_slim <- read_rds(here("Saved Data", 
                         "cleaned_length_at_age_raw.rds"))

glimpse(lt_slim)
lt_slim %>% 
  group_by(cr) %>% 
  summarise(
    n = n()
  )

lt_slim <- lt_slim %>% 
  mutate(
    wt_g = if_else(is.na(wt_g), true = 
                     round((10 ^ (-5.218126 +  3.038077 * log10(tl_mm))), 
                           digits = 0),
                   false = wt_g), 
  )

han <- read_csv(here("Data", 
                     "lake_charr_life_history_metric_muir_2021.csv")) %>%
  pivot_longer(cols = -c(Metric, n), 
               names_to = "percentile",
               values_to = "value") %>% 
  mutate(
    percentile = as.factor(str_remove(percentile, "%"))
  ) %>% 
  janitor::clean_names()

han
# ---- plot the data -----


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45, size = 3, aes(colour = sex)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45, size = 3, aes(colour = basin)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")




lt_sum <- lt_slim %>% 
  summarise(
    mean_tl = mean(tl_mm),
    sem_tl = sd(tl_mm) / sqrt(n()),
    mean_fl = mean(fl_mm), 
    sem_fl = sd(fl_mm) / sqrt(n()),
    mean_wt = mean(wt_g, ), 
    sem_wt = sd(wt_g) / sqrt(n())
  )

lt_sum  

lt_slim

lt_pap_han <- lt_slim %>% 
  select(age_est, tl_mm) %>% 
  rename(age = age_est, 
         tl = tl_mm) %>% 
  mutate(
    lake_name = "Papineau", 
    linf = NA, 
    t0 = NA, 
    k = NA
  ) %>% 
  select(lake_name, linf:k, age, tl) %>% 
  bind_rows(pred_vb, .)



# write.xlsx(lt_sum, here("Results", 
#                         "summary_morpho_age_fish.xlsx"))



# ---- compile data into a list for JAGS -----

pred_age <- seq(0, 26, 0.5) # adjust for plotting this is like tidy::crossing()

# data to be fed to JAGS as a list 
jags_data <- list(
  n_obs = nrow(lt_pap_han),
  length = lt_pap_han$tl,
  age = lt_pap_han$age,
  pred_age = pred_age,
  n_pred = length(pred_age)
)

# ---- Specify JAGS model code -----
# here you will create your priors based on known populations of lakers 
# from the literature to super von B. parameters 

jags_model <- function() {
  # priors
  
  linf ~ dunif(100, 800) 
  k ~ dunif(0, 0.5) 
  t0 ~ dunif(-1.25, 0.5)
  sig ~ dunif(0, 2) # using a log likelihood so log back transformed 
  sigma <- exp(sig / 2) # using a log likelihood so log back transformed 
  tau <- pow(sigma, -2)
  p ~ dunif(0, 1)
  # Likelihood 
  
  for (i in 1:n_obs) { 
    
    length[i] ~ dlnorm(log(length_hat[i]), 1 / sig ^ 2)
    
    length_hat[i] <- (linf) * (1 - exp(-k *(age[i] - t0))) 
    w[i] ~ dexp(tau)
    me[i] <- (1 - 2 * p) / (p * (1 - p)) * w[i] + length_hat[i]
    pe[i] <- (p * (1 - p) * tau) / (2 * w[i])
    y[i]  ~ dnorm(me[i], pe[i])
  }
  # Derived Quantities 
  # this makes it easier to plot 
  for (i in 1:n_pred){
    # expected length @ age 
    pred_length[i, 1] <- linf *  (1 - exp(-k *(pred_age[i] - t0))) 
    
    # sample random residuals to obtain random predictions 
    # when summarized this will represent the length @ age distr
    # that we can use to look at population density for length @ age 
    # 
    epi[i] ~ dlnorm(0, 1 / sig ^ 2)
    rand_length[i, 1] <- epi[i] * pred_length[i, 1]
    
  }
}

# ---- write model to a text file  to fed to JAGS -----
jags_file <- here("JAGS models", 
                  "lkt_von_b_bays_estimate_quantile.txt"
)


# use write_model which is from {postpack} and performs the same 
# opperation as R2OpenBUGS::write.model

write_model(jags_model, jags_file)

# ---- specify theta values for model -----
jags_inits <- function(nc) {
  inits = list()
  for (c in 1:nc) {
    inits[[c]] = list(
      linf = runif(1, 400, 800), # random value between 400 - 800 mm 
      k = runif(1, 0.05, 0.5), # random value between 0.1-0.5
      t0 = runif(1, -1.25, 0.5) # random value between 0-2 age-0
    )
  }
  return(inits)
}

# ---- identify and set which model parameters/nodes to monitor ----- 
jags_params <- c("linf", 
                 "k", 
                 "t0", 
                 "sig", 
                 "pred_length", 
                 "rand_length")

# ---- set MCMC dimensions -----
jags_dims = c(
  ni = 65000,  # number of post-burn-in samples per chain increased to 70000
  nb = 65000,  # number of burn-in samples increased to 70000
  nt = 20,     # thinning rate
  nc = 2 # number of chains
)

# ---- run the model with JAGS ----
post <- jagsUI::jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = jags_inits(jags_dims["nc"]),
  parameters.to.save = jags_params,
  n.adapt = 1000,
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = TRUE
)

# ---- model diagnostics with {postpack} ----- 

# nodes for diagnostics
diag_p <- c("linf", "k", "t0", "sig")

