# ---- load packages -----
library(brms)
library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(tidybayes)
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# ---- bring in data -----
lt_slim <- read_rds(here("Saved Data", 
                         "cleaned_length_at_age_raw.rds"))

glimpse(lt_slim)

lt_slim <- lt_slim %>% 
  rename(age = age_est) %>% 
  filter(age != is.na(age)) %>% 
  mutate(age = as.numeric(age))

# ---- set pirors -----
max(lt_slim$tl_mm)
mean(lt_slim$tl_mm)

prior_vis <- data.frame(
  param = rep(c('linf', 'k1', 't1', 'k2', 't2'), each = 1000),
  val = c(rnorm(1000, 760, 475),
          rnorm(1000, 0.119, 0.1),
          rnorm(1000, -0.783, 1),
          rnorm(1000, 0.097, 0.1),
          rnorm(1000, 1.01, 1))
)

ggplot(data = prior_vis, aes(x = val)) +
  geom_histogram(bins = 10) +
  facet_wrap(~param, scales = 'free_x')


vb_prior <- c(
  set_prior("normal(740, 475)", nlpar = 'Linf', lb = 700),
  set_prior("normal(0.12, 0.05)", nlpar = 'k1', lb = 0),
  set_prior("normal(0.099, 0.05)", nlpar = 'k2', lb = 0),
  set_prior("normal(-0.583, 1)", nlpar = 't1'),
  set_prior("normal(-2.01, 1)", nlpar = 't2')
)

vonb <- bf(
  # formula of the equation we want to fit
  tl_mm ~ (age < (k2*t2 - k1*t1) / (k2 - k1)) * (Linf * (1 - exp(-k1 * (age - t1)))) +
    (age >= (k2*t2 - k1*t1) / (k2 - k1)) * (Linf * (1 - exp(-k2 * (age - t2)))),
  
  # Parameters we want to fit
  #   In the future we can create "random effects" by grouping the parameters
  #   Linf + k2 + t2 ~ 1 | lab_group for example
  Linf + k1 + k2 + t1 + t2 ~ 1,
  
  # Tell brms that this is a non-linear equation
  nl = T
)



vbfit <- brm(vonb, data = lt_slim, prior = vb_prior,
             iter = 5000, cores = 4, thin = 5)

beepr::beep()
plot(vbfit, N = 2, ask = F)


pp_check(vbfit)
pp_check(vbfit, type = 'ecdf_overlay')
summary(vbfit)

plot(conditional_effects(vbfit), points = T)
tp <- vbfit |>
  spread_draws(`b_.*`, regex = TRUE) |> 
  rename_all(function(.) gsub('[b_]|Intercept', '', .)) |> 
  mutate(tp = (k2*t2 - k1*t1) / (k2 - k1))

quantile(tp$tp, c(0.025, 0.5, 0.95))
         