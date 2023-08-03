# ---- bring in packages -----
{
  library(broom)
  library(broom.mixed)
  library(car)
  library(data.table)
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(here)
  library(lubridate)
  library(lemon)
  library(lme4)
  library(mgcv)
  library(quantreg)
  library(openxlsx)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
}


# ---- bring data ---- 
han_raw <- read_csv(here("Data",  
                         "hansen_2021_raw_data.csv")) %>% 
  janitor::clean_names()

# ---- select vb curve parameters and remove rows with NA -----
vb <- han_raw %>% 
  select(lake_name, linf, t0, k)

vb_all <- vb 

sum_n <- vb_all %>% 
  pivot_longer(cols = -lake_name, 
               names_to = "pram", 
               values_to = "est") %>% 
  group_by(pram) %>% 
  summarise(
    n = sum(!is.na(est))
  ) %>% 
  ungroup()
sum_n
# ---- bring papineau predicted length @ age and VB coeff -----
vb_pl_wide <- read_rds(here("Saved Data", 
                            "bayes_predicted_model_wide.rds"))

vb_pl_coef <- read_rds(here("Saved data", 
                            "von_b_50_coef.rds"))
# ---- modeify VB coef ----
vb_pl_coef <- vb_pl_coef %>% 
  mutate(
    lake_name = "Papineau"
  ) %>% 
  select(lake_name, linf, t0, k)

# ---- join all for quantile reggression ----- 
vb_all <- bind_rows(vb_all, 
                    vb_pl_coef)
# add in row id that is lake id 
vb_all <- vb_all %>% 
  mutate(
    id = 1:nrow(.), 
    lake_name = as.factor(lake_name)
  )


ggplot(data = vb_all, aes(x = id, y = linf)) + 
  geom_point(size = 4) + 
  geom_quantile(quantiles = seq(0.25, 0.75, 0.25), show.legend = TRUE)

ggplot(data = vb_all, aes(x = linf)) + 
  geom_histogram(bins = 14, fill = "white", colour = "black") + 
  # geom_vline(xintercept = c(
  #   linf_qr$coefficients[1],
  #   linf_qr$coefficients[3],
  #   linf_qr$coefficients[5]), linetype = 1, 
  #   colour = "blue", 
  #   linewidth = 1
  # )+ 
  # geom_vline(xintercept = c(
  #   605,
  #   717,
  #   880), linetype = 4, 
  #   colour = "red", 
  #   linewidth = 1
  # )+ 
  scale_x_continuous(breaks = seq(400, 1600, 200)) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank()
  )

# linf = c(
#   linf_qr$coefficients[1],
#   linf_qr$coefficients[3],
#   linf_qr$coefficients[5]
# 
# vb_all$random <- runif(489, min = 1, max = 489)
# # ---- use quantile regressions for each coeff ---- 
# linf_qr <- rq(linf ~ random, data = vb_all,
#               tau = seq(0.25, 0.75, 0.25))
# linf_qr
# 
# summary(linf_qr)
# 
# 
# 
# k_qr <- rq(k ~ id, data = vb_all, tau = seq(0.25, 0.75, 0.25))
# t0_qr <- rq(t0 ~ id, data = vb_all, tau = seq(0.25, 0.75, 0.25))
# 
# rq()

quantile(vb_all$linf, na.rm = TRUE)
quantile(vb_all$k, na.rm = TRUE)
quantile(vb_all$t0, na.rm = TRUE)

# qr_vb <- tibble(
#   qr = factor(seq(25, 75, 25)), 
#   linf = c(
#     linf_qr$coefficients[1],
#     linf_qr$coefficients[3],
#     linf_qr$coefficients[5]
#   ),
#   t0 = c(
#     t0_qr$coefficients[1],
#     t0_qr$coefficients[3],
#     t0_qr$coefficients[5]
#   ),
#   k = c(
#     k_qr$coefficients[1],
#     k_qr$coefficients[3],
#     k_qr$coefficients[5]
#   )
# ) %>% 
#   expand_grid(
#     .,
#     age = seq(0, 26, 0.5)
#   ) %>% 
#   mutate(
#     # age_est = (log(1 - (tl / l_inf)) / -k ) + t0
#     tl = linf * (1 - exp(-k *(age - t0))) 
#   )



ggplot() + 
  geom_line(data = qr_vb, aes(x = age, y = tl, linetype = qr, group = qr)) + 
  geom_line(data = vb_pl_wide, aes(x = pred_age, y = m_50), linewidth = 1) + 
  scale_linetype_manual(
    name = "Percentile", 
    values = c(2:4)
  ) 
summary(linf_qr)

summary(k_qr)
summary(linf_qr)
summary(t0_qr)

vb_all_p <- expand_grid(
  vb_all,
  age = seq(0, 26, 0.5)
) %>% 
  mutate(
    # age_est = (log(1 - (tl / l_inf)) / -k ) + t0
    tl = linf * (1 - exp(-k *(age - t0))) 
  )
pred_vb
# ---- preliineary plot 
pred_vb <- expand_grid(
  vb_all,
  age = seq(0, 26, 0.5)
) %>% 
  mutate(
    # age_est = (log(1 - (tl / l_inf)) / -k ) + t0
    tl = linf * (1 - exp(-k *(age - t0))) 
  )
ggplot() + 
  geom_line(data = pred_vb, aes(x = age, y = tl, 
                                group = lake_name,
                                # colour = lake_name
  )) + 
  geom_line(data = vb_pl_wide, aes(x = pred_age, y = m_50), 
            colour = "black", linewidth = 1, linetype = 2) + 
  scale_colour_viridis_d(end = 0.85, name = "Lakes")





pred_vb
