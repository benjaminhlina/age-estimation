# ---- load packages -----

library(data.table)
library(dplyr)
library(here)
library(forcats)
library(ggplot2)
library(ggridges)
library(ggExtra)
library(ggmcmc)
library(postpack)
library(patchwork)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

# ---- read data -----

lt_slim <- read_rds(here("Saved Data", 
                         "cleaned_length_at_age_raw.rds"))

glimpse(lt_slim)
# ---- plot the data -----


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45, size = 3, aes(colour = sex)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")


lt_slim <- lt_slim %>% 
  filter(basin != is.na(basin)) %>% 
  mutate(
    basin = factor(stringr::str_replace(as.character(basin), " Basin",
                                        ""),
                   levels = c("East", "West", "North")),
    basin_n = as.numeric(basin)
  )


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45, size = 3, aes(colour = basin_n)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")



lt_slim %>% 
  group_by(basin) %>% 
  summarise(n = n ()) %>% 
  ungroup()

# ---- compile data into a list for JAGS -----

pred_age <- seq(0, 26, 0.5) # adjust for plotting this is like tidy::crossing()

# data to be fed to JAGS as a list 
jags_data <- list(
  n_obs = nrow(lt_slim),
  basin = lt_slim$basin_n, 
  length = lt_slim$tl_mm,
  age = lt_slim$age_est,
  pred_age = pred_age,
  n_pred = length(pred_age)
)

jags_data$n_obs

# ---- Specify JAGS model code -----
# here you will create your priors based on known populations of lakers 
# from the literature to super von B. parameters 

jags_model <- function() {
  # priors
  linf ~ dunif(100, 800) 
  k ~ dunif(0, 0.5) 
  t0 ~ dunif(-1.25, 0.5)
  sig ~ dunif(0, 2) # using a log likelihood so log back transformed 
  linf_diff ~ dunif(-500, 500) # estimated quantity for lenght @ infiity 
  k_diff ~ dunif(-0.05, 0.05) 
  t0_diff ~ dunif(-0.2, 0.2)
  # Likelihood 
  
  for (i in 1:n_obs) { 
    
    length[i] ~ dlnorm(log(length_hat[i]), (1 / sig ^ 2))
    length_2[i] ~ dlnorm(log(length_hat_2[i]), (1 / sig ^ 2))
    length_3[i] ~ dlnorm(log(length_hat_3[i]), (1 / sig ^ 2))
    # length[i] ~ dlnorm(log(length_hat[i]), (1 / sig ^ 2))
    # length[i] ~ dlnorm(log(length_hat[i]), 1 / sig ^ 2)
    
    # length_hat[i] <- (linf) * (1 - exp(-k *(age[i] - t0))) 
    # length_hat[i] <- (linf + linf_diff * basin[i]) * (1 - exp(-k *(age[i] - t0))) 
    length_hat[i] <- (linf + linf_diff * basin[i]) * (1 - exp((-k * (age[i] - t0)))) 
    length_hat_2[i] <- (linf) * (1 - exp((-k + k_diff * basin[i]) * (age[i] - t0 ))) 
    length_hat_3[i] <- (linf) * (1 - exp(-k * (age[i] - (t0 + t0_diff * basin [i])))) 
  }
  # Derived Quantities 
  # this makes it easier to plot 
  for (i in 1:n_pred){
    # expected length @ age 
    pred_length[i, 1] <- linf *  (1 - exp(-k *(pred_age[i] - t0))) 
    pred_length[i, 2] <- (linf + linf_diff) * 
      (1 - exp((-k + k_diff) *(pred_age[i] - (t0 + t0_diff))))
    pred_length[i, 3] <- (linf - linf_diff) * 
      (1 - exp((-k - k_diff) *(pred_age[i] - (t0 - t0_diff))))
    # sample random residuals to obtain random predictions 
    # when summarized this will represent the length @ age distr
    # that we can use to look at population density for length @ age 
    # 
    epi[i] ~ dlnorm(0, 1 / sig ^ 2)
    rand_length[i, 1] <- epi[i] * pred_length[i, 1]
    rand_length[i, 2] <- epi[i] * pred_length[i, 2]
    rand_length[i, 3] <- epi[i] * pred_length[i, 3]
  }
}

# ---- write model to a text file  to fed to JAGS -----
jags_file <- here("JAGS models", 
                  "lkt_von_b_bays_estimate_basin.txt"
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
      linf_diff = runif(1, -50, 50), # random value between -50 and 50 to start
      # differences 
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
                 "linf_diff",
                 "k_diff",
                 "t0_diff",
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
diag_p <- c("linf", "k", "t0", "sig", "linf_diff", "k_diff", "t0_diff")

# view convergence diagnostic summaries
post_summ(post, diag_p, Rhat = TRUE, neff = TRUE)[c("Rhat", "neff"),]

# view diagnostic plots
gg_post <- ggs(post) %>%
  janitor::clean_names() %>%
  filter(parameter %in% c("k", "linf",
                          "sig", "t0", "linf_diff", 
                          "k_diff", "t0_diff")) 
# gg_post_all <- ggs(post)

ggplot(data = gg_post, aes(x = value)) + 
  geom_density(aes(fill = as.factor(chain), 
                   colour = as.factor(chain)), 
               linewidth = 1.2) +
  scale_fill_viridis_d(alpha = 0.25,
                       name = "Chain",
                       end = 0.85, begin = 0.35, option = "D") +
  scale_colour_viridis_d(alpha = 0.25,
                         name = "Chain",
                         end = 0.65, begin = 0.35, option = "D") + 
  facet_wrap(~ parameter, ncol = 1, scales = "free") + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.955, 0.945), 
  ) + 
  labs(x = "Estimate", 
       y = "Density") + 
  
  
  ggplot(data = gg_post, aes(x = iteration, y = value)) + 
  geom_line(aes(colour = as.factor(chain)), 
  ) +
  scale_fill_viridis_d(
    alpha = 0.7,
    name = "Chain",
    end = 0.85, begin = 0.35, option = "D") +
  scale_colour_viridis_d(
    alpha = 0.7,
    name = "Chain",
    end = 0.65, begin = 0.35, option = "D") + 
  facet_wrap(~ parameter, ncol = 1, scales = "free") + 
  scale_x_continuous(labels = seq(attributes(gg_post)$nBurnin,
                                  125000,
                                  length.out = 4)) +
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.955, 0.945), 
    legend.background = element_blank()
  ) + 
  labs(x = "Iteration", 
       y = "Estimate") -> p_diag

p_diag

# diag_plots(post, diag_p, ext_device = TRUE)


gg_posts <- ggs(post) %>%
  # janitor::clean_names() %>%
  filter(Parameter %in% c("k", "linf",
                          "sig", "t0", "linf_diff", 
                          "k_diff", "t0_diff"))
# ggs_traceplot(gg_posts)
ggs_running(gg_posts)

ggs_compare_partial(gg_posts)
ggs_autocorrelation(gg_posts) 

# ---- what's the probability the difference in linf ---- 
# is it positive or negative?
linf_diff <- as.numeric(post_subset(post, "linf_diff", TRUE))
p_neg <- mean(linf_diff < 0)
p_neg # posterior probablity of east basin vs west basin is 
p_pos <- mean(linf_diff > 0)
p_pos


# what's the probability the difference in linf is negative, but "negligible"?
mean(linf_diff > -10 & linf_diff < 0)
# ---- look at model parameter estimates -----
post_summ(post, diag_p)

post_summ_coef <- cbind(metric = post_summ(post, diag_p) %>% 
                          row.names(.), 
                        post_summ(post, diag_p) %>% 
                          as_tibble(.)) %>% 
  tibble()
# post_summ_coef_25 <- post_summ_coef
# post_summ_coef_25
post_summ_coef
# post_summ_coef <- post_summ_coef %>% 
#   mutate(
#     w = linf *k
#   )
# post_summ_coef
# write_rds(post_summ_coef, here("Saved Data",
#                                "von_b_bays_coef.rds"))


# ---- create and manipulate predicted length for plotting -----

pred_length_1 <- post_summ(post, "pred_length[.+,1]") # regexpres
pred_length_2 <- post_summ(post, "pred_length[.+,2]") # regexpres
pred_length_3 <- post_summ(post, "pred_length[.+,3]") # regexpres



# tke matrix and create tibble 
pred_length_df <- bind_rows(.id = "basin", 
                            do.call(bind_rows, 
                                    lapply(pred_length_1, function(x) as_tibble(x))),
                            do.call(bind_rows, 
                                    lapply(pred_length_2, function(x) as_tibble(x))),
                            do.call(bind_rows, 
                                    lapply(pred_length_3, function(x) as_tibble(x))),
)
pred_length_df <- cbind(metric = row.names(pred_length_1), 
                        pred_length_df) %>% 
  tibble() %>% 
  mutate(
    metric = factor(case_when(
      metric %in% "mean" ~ "mean", 
      metric %in% "sd" ~ "sd", 
      metric %in% "50%" ~ "m_50",
      metric %in% "2.5%" ~ "m_2.5", 
      metric %in% "97.5%" ~ "m_97.5"), 
      level = c("mean",
                "sd",
                "m_50",
                "m_2.5", "m_97.5")), 
    basin = factor(case_when(
      basin == 1 ~ "East", 
      basin == 2 ~ "West", 
      basin == 3 ~ "North"),
      levels = c("East", "West", "North")),
  )

# makee tibble wide 
pred_lenght_df_wide <- pred_length_df %>% 
  group_by(metric, basin) %>% 
  mutate(row = row_number(), 
         pred_age =  seq(0, 26, 0.5)) %>%
  pivot_wider(names_from = "metric", 
              values_from = "value") 

# grab just the model parameters used to predict 
post_summ_coef %>% 
  filter(metric %in% "50%") -> model_fit

mfl <- model_fit %>% 
  pivot_longer(cols = -metric, 
               names_to = "coefficient",
               values_to = "est")
mfl
model_fit %>% 
  write_rds(here("Saved data", 
                 "von_b_50_coef.rds"))


# extract each parmater indivually to put von B. curve exquation on plot

linf <- formatC(model_fit$linf, format="f", digits = 0)
k <- formatC(model_fit$k, format = "f", digits = 3)
# Handle t0 differently because of minus in the equation
t0 <- model_fit$t0
t0 <- paste0(ifelse(t0 < 0, "+", "-"), formatC(abs(t0), format = "f", 
                                               digits = 3))

linf_diff_l <- formatC(model_fit$linf + model_fit$linf_diff , format="f", digits = 0)
k_diff <- formatC(model_fit$k + model_fit$k_diff, format = "f", digits = 3) 
t0_diff <- model_fit$t0 + model_fit$t0_diff
t0_diff <- paste0(ifelse((t0_diff) < 0, "+", "-"), formatC(abs(t0_diff), format = "f", 
                                                           digits = 3))

linf_diff_n <- formatC(model_fit$linf - model_fit$linf_diff , format="f", digits = 0)
k_diff_n <- formatC(model_fit$k - model_fit$k_diff, format = "f", digits = 3) 
t0_diff_n <- model_fit$t0 - model_fit$t0_diff
t0_diff_n <- paste0(ifelse((t0_diff_n) < 0, "+", "-"), formatC(abs(t0_diff_n), format = "f", 
                                                             digits = 3))
# Put together and return
# labels <- c(
#   paste0("TL ==", linf_diff_l,"~bgroup('(',1-e^{-", k_diff,"~(age", t0_diff,")},')')"), 
#   paste0("TL ==", linf,"~bgroup('(',1-e^{-", k,"~(age", t0,")},')')"),
#   paste0("TL ==", linf_diff_n,"~bgroup('(',1-e^{-", k_diff_n,"~(age", t0_diff_n,")},')')") 
# )
label_1 <- c(paste0("TL ==", linf_diff_l,"~bgroup('(',1-e^{-", 
                    k_diff,"~(age", t0_diff,")},')')"))

label_2<- c(paste0("TL ==", linf,"~bgroup('(',1-e^{-", k,"~(age", t0,")},')')"))
label_3 <- c(paste0("TL ==", linf_diff_n,"~bgroup('(',1-e^{-", 
                    k_diff_n,"~(age", t0_diff_n,")},')')")) 


# to grab the correct colour 
# alpha_col <- function(color_name, alpha = 0.35) {
#   rgb <- col2rgb(color_name)
#   rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = (1 - alpha) * 255)
# }
# ---- plot Bayes Von B curve using ggplot -----  
ggplot() + 
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm, colour = basin),
             alpha = 0.45, size = 3) +
  geom_line(data = pred_lenght_df_wide, 
            aes(x = pred_age, y = m_50, colour = basin), 
            linewidth = 1) + 
  geom_ribbon(data = pred_lenght_df_wide,
              aes(ymin = m_2.5,
                  ymax = m_97.5,
                  x = pred_age, y = m_50, fill = basin), 
              alpha = 0.15) +  
  scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", name = "Basin") + 
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "Basin") + 
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(0, 700, 100)) + 
  annotate(geom = "text",
           label = label_2,
           parse = TRUE,
           size = 4,
           x = 22.75,  
           y = 85, 
           hjust = 1.1,
           vjust = -0.5
  ) +
  annotate(geom = "text",
           label = label_1,
           parse = TRUE,
           size = 4,
           x = 22.75, 
           y = 25, 
           hjust = 1.1,
           vjust = -0.5
  ) +
  annotate(geom = "text",
           label = label_3,
           parse = TRUE,
           size = 4,
           x = 22.75, 
           y = -30, 
           hjust = 1.1,
           vjust = -0.5
  ) +
  theme_classic(base_size = 15) + 
  theme(
    legend.position = c(0.075, 0.90), 
    # legend.background = element_blank()
  ) + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)") -> p



p

# ggsave(filename = here("Plots",
#                        "von B Bayes",
#                        "growth_curve_lt_bayes_basin_diff.png"), 
#        plot = p, width = 11, height = 8.5)

# p
write_rds(p, here("Saved Plots", 
                  "growth_curve_lt_bayes_basin.rds"))
# hypothetically speaking, if a fish matures at 350 mm in papineau what age 

# l_crit <- 350
# par(mar = c(3, 3, 1, 1), tcl = -0.25, mgp = c(2, 0.5, 0))
# rand_length <- post_subset(post, "rand_length[.+,1]", TRUE)
# p_crit <- apply(rand_length, 2, function(x) mean(x > l_crit))
# 
# plot(p_crit ~ pred_age, 
#      ylab = "Pr(Mature)", 
#      xlab = "Age (Years)", 
#      type = "l", col = "blue", ylim = c(0, 1))


# ---- plot prep for density distribution of ages -----

# sample random residuals to obtain random predictions 
# when summarized this will represent the length @ age distr
# that we can use to look at population density for length @ age 
# 
# select every 2 year increments 
keep_pred_age <- which(pred_age %in% seq(0, 26, 2))
# extract random lengths from predicted model 
rand_length_1 <- post_subset(post, "rand_length[.+,1]", 
                             matrix = TRUE
)
rand_length_2 <- post_subset(post, "rand_length[.+,2]", 
                             matrix = TRUE
)
rand_length_3 <- post_subset(post, "rand_length[.+,3]", 
                             matrix = TRUE
)



# create random length per 2 yr age category for density ploting
rand_length_df <- bind_rows(.id = "basin",
                            as_tibble(rand_length_1[,keep_pred_age]) %>% 
                              mutate(id = 1:nrow(.)), 
                            as_tibble(rand_length_2[,keep_pred_age]) %>% 
                              mutate(id = 1:nrow(.)), 
                            as_tibble(rand_length_3[,keep_pred_age]) %>% 
                              mutate(id = 1:nrow(.)), 
) %>% 
  pivot_longer(cols = -c("id", "basin"), 
               names_to = "age", 
               values_to = "rand_length") %>%  
  mutate(
    age_pred = case_when(
      age %in% "rand_length[1,1]" ~ 0,
      age %in% "rand_length[5,1]" ~ 2, 
      age %in% "rand_length[9,1]" ~ 4, 
      age %in% "rand_length[13,1]" ~ 6, 
      age %in% "rand_length[17,1]" ~ 8, 
      age %in% "rand_length[21,1]" ~ 10, 
      age %in% "rand_length[25,1]" ~ 12, 
      age %in% "rand_length[29,1]" ~ 14, 
      age %in% "rand_length[33,1]" ~ 16, 
      age %in% "rand_length[37,1]" ~ 18, 
      age %in% "rand_length[41,1]" ~ 20, 
      age %in% "rand_length[45,1]" ~ 22, 
      age %in% "rand_length[49,1]" ~ 24, 
      age %in% "rand_length[53,1]" ~ 26, 
      age %in% "rand_length[1,2]" ~ 0,
      age %in% "rand_length[5,2]" ~ 2, 
      age %in% "rand_length[9,2]" ~ 4, 
      age %in% "rand_length[13,2]" ~ 6, 
      age %in% "rand_length[17,2]" ~ 8, 
      age %in% "rand_length[21,2]" ~ 10, 
      age %in% "rand_length[25,2]" ~ 12, 
      age %in% "rand_length[29,2]" ~ 14, 
      age %in% "rand_length[33,2]" ~ 16, 
      age %in% "rand_length[37,2]" ~ 18, 
      age %in% "rand_length[41,2]" ~ 20, 
      age %in% "rand_length[45,2]" ~ 22, 
      age %in% "rand_length[49,2]" ~ 24, 
      age %in% "rand_length[53,2]" ~ 26, 
      age %in% "rand_length[1,3]" ~ 0,
      age %in% "rand_length[5,3]" ~ 2, 
      age %in% "rand_length[9,3]" ~ 4, 
      age %in% "rand_length[13,3]" ~ 6, 
      age %in% "rand_length[17,3]" ~ 8, 
      age %in% "rand_length[21,3]" ~ 10, 
      age %in% "rand_length[25,3]" ~ 12, 
      age %in% "rand_length[29,3]" ~ 14, 
      age %in% "rand_length[33,3]" ~ 16, 
      age %in% "rand_length[37,3]" ~ 18, 
      age %in% "rand_length[41,3]" ~ 20, 
      age %in% "rand_length[45,3]" ~ 22, 
      age %in% "rand_length[49,3]" ~ 24, 
      age %in% "rand_length[53,3]" ~ 26, 
    ), 
    basin = factor(case_when(
      basin == 1 ~ "East", 
      basin == 2 ~ "West", 
      basin == 3 ~ "North"),
      levels = c("East", "West", "North")), 
  ) %>% 
  select(basin, age_pred, rand_length) %>% 
  arrange(age_pred) %>% 
  mutate(age_pred_f = as.factor(age_pred)) %>% 
  drop_na()

rand_length_df



# ---- plot densities -----
ggplot(rand_length_df, aes(x = rand_length)) + 
  geom_density(aes(
    colour = age_pred_f
    # fill = age_pred_f
  ), alpha = 0.25, 
  linewidth = 0.8, 
  # colour = "black"
  ) + 
  facet_wrap(~ basin, ncol = 1) + 
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", "Estimated Age (Yr)") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85,
                       option = "D", "Estimated Age (Yr)") +
  theme_bw(base_size = 15) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) + 
  theme(panel.grid = element_blank(), 
        # legend.position = c(0.88, 0.86),
        legend.background = element_blank(), 
        axis.text = element_text(colour = "black"), 
  ) +
  guides(
    fill = guide_legend(ncol = 2),
    colour = guide_legend(ncol = 2)
  ) + 
  labs(x = "Total Length (mm)", 
       y = "Density") -> p3
p3
ggsave(filename = here("Plots", 
                       "posterior distribution",
                       "predicted_density_lkt_age_length_basin_fill_diff.png"), 
       plot = p3, height = 8.5, width = 11)
write_rds(p3, here("Saved Plots", 
                  "posterior_distribution_lkt_age_length_basin_.rds"))
# facet densities 
ggplot(rand_length_df, aes(x = rand_length)) + 
  geom_density(
    # aes(colour = age_pred_f), 
    linewidth = 1) + 
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", "Estimated Age (Yr)") + 
  theme_bw(base_size = 15) + 
  facet_wrap(. ~ age_pred_f) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) +
  
  # theme(panel.grid = element_blank(), 
  #       legend.position = c(0.92, 0.9)) +
  guides(colour = guide_legend(ncol = 4)) + 
  labs(x = "Total Length (mm)", 
       y = "Density")

# rand_length_df  <- rand_length_df %>% 
#   mutate(
#     age_pred_f = fct_rev(age_pred_f)
#   )
# ---- Tornado plot of posterior densities for random residuals ---- 
# posterior distribution of length-at-age 
# (including between-individual variability)
ggplot(rand_length_df, aes(x = rand_length)) + 
  geom_density_ridges(
    aes(
      y = age_pred_f, 
      fill = age_pred_f), 
    # linewidth = 1
  ) + 
  facet_wrap(~ basin, ncol = 1) + 
  scale_x_continuous(breaks = seq(0, 1200, 100)) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, alpha = 0.75,
                       option = "D", "Estimated Age (Yr)") + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none") +
  guides(fill = guide_legend(ncol = 4)) + 
  labs(x = "Total Length (mm)", 
       y = "Estimated Age (Yr)") -> p4
p4
# "predicted_density_ridge_lkt_age_length"
ggsave(filename = here("Plots",
                       "posterior distribution",
                       "predicted_density_ridge_tornado_age_length_basin_diff.png"), 
       plot = p4, height = 8.5, width = 11)

write_rds(p4, here("Saved Plots", 
                   "redicted_density_ridge_tornado_age_length_basin_diff.rds"))

# ---- prep to estimate age from a given length -----
# 
# 
# 
# vcov_mcmc <- dclone::vcov.mcmc.list(post)
# vcov_matrix <- vcov_mcmc[-4:-109, -4:-109]
# vcov_matrix <- vcov_matrix[-1:-1, -1:-1]
# 
# mfl <- mfl %>% 
#   mutate(
#     coefficient = factor(coefficient, 
#                          level = c("k", "linf", "sig", "t0"))
#   ) %>% 
#   arrange(coefficient)
# mfl
# model_fit
# 
# 
# hill <- MASS::mvrnorm(n = 10000, mu = mfl$est, Sigma = vcov_matrix)
# hill <- data.frame(hill)
# 
# head(hill)
# 
# 
# { ggplot(data = hill, aes(x = t0, y = k)) +
#     geom_density2d_filled(show.legend = F)+
#     geom_point(alpha = 0) +
#     annotate('point', x = model_fit$t0, y = model_fit$k)
# } |> 
#   ggMarginal()
# 
# 
# { ggplot(data = hill, aes(x = linf, y = k)) +
#     geom_density2d_filled(show.legend = F)+
#     geom_point(alpha = 0) +
#     annotate('point', x = model_fit$linf, y = model_fit$k)} |> 
#   ggMarginal()
# 
# 
# { ggplot(data = hill, aes(x = linf, y = t0)) +
#     geom_density2d_filled(show.legend = F)+
#     geom_point(alpha = 0) +
#     annotate('point', x = model_fit$linf, y = model_fit$t0)} |> 
#   ggMarginal()
# 
# 
# pred_age <- function(tl, linf = 2241, k = 0.097, t0 = 1.09) {
#   (log(1 - (tl / linf)) / -k ) + t0
# }
# 
# hill <- hill |>
#   mutate(pred_age = pred_age(tl =  485, t0 = t0, linf = linf, k = k)) 
# glimpse(hill)
# 
# ggplot(data = hill) +
#   geom_density(aes(x = pred_age)) +
#   geom_vline(xintercept = 11.677702)+
#   geom_vline(xintercept = quantile(hill$pred_age, c(0.025, 0.975), na.rm = T), color = 'red') +
#   xlim(6, quantile(hill$pred_age, 0.99, na.rm = T))
# 
# 
# pred_age_pp <- function(., quant) {suppressWarnings(
#   quantile(pred_age(., t0 = hill$t0,
#                     linf = hill$linf,
#                     k = hill$k), probs = quant, na.rm = T)
# )}
# 
# 
# df <- read_rds(here("Saved data", 
#                     "tag_iso_est_age.rds"))
# 
# df
# df$pred_age_med <- sapply(df$tl_mm, pred_age_pp, quant = 0.5)
# df$lpi_age <- sapply(df$tl_mm, pred_age_pp, quant = 0.025)
# df$upi_age <- sapply(df$tl_mm, pred_age_pp, quant = 0.975)
# glimpse(df)
# 
# 
# ggplot(data = df) +
#   geom_pointrange(aes(x = age_est, y = tl_mm, 
#                       xmin = lpi_age, xmax = upi_age,
#                       group = sample)) +
#   scale_x_continuous(breaks = seq(0, 80, 5)) + 
#   geom_point(aes(x = pred_age_med, y = tl_mm,
#                  group = sample), shape = 'triangle') +
#   theme_bw(base_size = 15) + 
#   
#   labs(x = 'Predicted age', 
#        y = "Total Length (mm)") -> p1
# 
# p1
# ggsave(filename = here("Plots",
#                        "predicted age for length SI",
#                        "predicted_age_for_length_SI.png"), 
#        plot = p1, height = 8.5, width = 11)
