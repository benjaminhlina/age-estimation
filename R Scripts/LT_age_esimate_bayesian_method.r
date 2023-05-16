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


  ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45, size = 3, aes(colour = basin)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")

# ---- compile data into a list for JAGS -----

pred_age <- seq(0, 26, 0.5) # adjust for plotting this is like tidy::crossing()

# data to be fed to JAGS as a list 
jags_data <- list(
  n_obs = nrow(lt_slim),
  length = lt_slim$tl_mm,
  age = lt_slim$age_est,
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
  
  # Likelihood 
  
  for (i in 1:n_obs) { 
    
    length[i] ~ dlnorm(log(length_hat[i]), 1 / sig ^ 2)
    
    length_hat[i] <- (linf) * (1 - exp(-k *(age[i] - t0))) 
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
                  "lkt_von_b_bays_estimate.txt"
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

# view convergence diagnostic summaries
post_summ(post, diag_p, Rhat = TRUE, neff = TRUE)[c("Rhat", "neff"),]

# view diagnostic plots
gg_post <- ggs(post) %>%
  janitor::clean_names() %>%
  filter(parameter %in% c("k", "linf",
                          "sig", "t0")) 
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
                          "sig", "t0"))
# ggs_traceplot(gg_posts)
ggs_running(gg_posts)

ggs_compare_partial(gg_posts)
ggs_autocorrelation(gg_posts) 
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

pred_length <- post_summ(post, "pred_length[.+,1]") # regexpres


# tke matrix and create tibble 
pred_length_df <- do.call(bind_rows, 
                          lapply(pred_length, function(x) as_tibble(x)))
pred_length_df <- cbind(metric = row.names(pred_length), 
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
                "m_2.5", "m_97.5")
    )
  )

# makee tibble wide 
pred_lenght_df_wide <- pred_length_df %>% 
  group_by(metric) %>% 
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
# extract each parmater indivually to put von B. curve exquation on plot
linf <- formatC(model_fit$linf, format="f", digits = 0)
k <- formatC(model_fit$k, format = "f", digits = 3)
# Handle t0 differently because of minus in the equation
t0 <- model_fit$t0
t0 <- paste0(ifelse(t0 < 0, "+", "-"), formatC(abs(t0), format = "f", 
                                               digits = 3))
# Put together and return
labels <- paste("TL ==", linf,"~bgroup('(',1-e^{-", k,"~(age", t0,")},')')")

# to grab the correct colour 
alpha_col <- function(color_name, alpha = 0.35) {
  rgb <- col2rgb(color_name)
  rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = (1 - alpha) * 255)
}
# ---- plot Bayes Von B curve using ggplot -----  
ggplot() + 
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm),
             alpha = 0.45, size = 3) +
  geom_line(data = pred_lenght_df_wide, 
            aes(x = pred_age, y = m_50), linewidth = 1, colour = "#39808f") + 
  geom_ribbon(data = pred_lenght_df_wide,
              aes(ymin = m_2.5,
                  ymax = m_97.5,
                  x = pred_age, y = m_50),
              fill = "#39808f", alpha = 0.25) +  
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(0, 700, 100)) + 
  annotate(geom = "text",
           label = labels,
           parse = TRUE,
           size = 4,
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -0.5) + 
  theme_classic(base_size = 15) + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)") -> p



p

ggsave(filename = here("Plots",
                       "von B Bayes",
                       "growth_curve_lt_bayes.png"), 
       plot = p, width = 11, height = 8.5)


write_rds(p, here("Saved Plots", 
                  "growth_curve_lt_bayes.rds"))
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
rand_length <- post_subset(post, "rand_length[.+,1]", 
                           matrix = TRUE
)

# create random length per 2 yr age category for density ploting
rand_length_df <- as_tibble(rand_length[,keep_pred_age]) %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(cols = -id, 
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
    )
  ) %>% 
  select(age_pred, rand_length) %>% 
  arrange(age_pred) %>% 
  mutate(age_pred_f = as.factor(age_pred))

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
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", "Estimated Age (Yr)") + 
  # scale_fill_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", "Estimated Age (Yr)") + 
  theme_bw(base_size = 15) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.88, 0.86),
        legend.background = element_blank(), 
        axis.text = element_text(colour = "black"), 
  ) +
  guides(
    # fill = guide_legend(ncol = 2)
    colour = guide_legend(ncol = 2)
    ) + 
  labs(x = "Total Length (mm)", 
       y = "Density") -> p3
p3
ggsave(filename = here("Plots", 
                       "posterior distribution",
                       "predicted_density_lkt_age_length.png"), 
       plot = p3, height = 8.5, width = 11)

write_rds(p3, here("Saved Plots", 
                  "posterior distribution_lkt_age_length.rds"))

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
                       "predicted_density_ridge_tornado_age_length.png"), 
       plot = p4, height = 8.5, width = 11)


write_rds(p4, here("Saved Plots", 
                   "predicted_density_ridge_tornado_age_length.rds"))
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
