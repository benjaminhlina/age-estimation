

# ---- load packages -----

library(data.table)
library(DHARMa)
library(dplyr)
library(here)
library(forcats)
library(ggplot2)
library(gratia)
library(emmeans)
library(ggridges)
library(ggExtra)
library(ggmcmc)
library(glmmTMB)
library(mgcv)
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

model_fit <-  read_rds(here("Saved Data", 
                            "vb_param_bayes.rds"))

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

pred_length_wide_df <- read_rds(here("Saved data", 
                                     "bayes_predicted_model_wide.rds"))
pred_length <- read_rds(here("Saved data", 
                             "bayes_predicted_model.rds"))
# ---- plot the data 

han_vb <- han %>% 
  filter(metric %in% c("L1", "t0", "K")) %>% 
  dplyr::select(-n)

han_vb_wide <- han_vb %>% 
  pivot_wider(
    names_from = metric, values_from = value
  ) %>% 
  rename(
    l_inf = L1, 
    k = K
  ) %>% 
  mutate(
    t0 = (as.numeric(gsub("\\–", "", t0)) * -1),
    k = as.numeric(k),
    l_inf = as.numeric(l_inf)
  )

han_vb_25_75 <- han_vb_wide %>% 
  filter(percentile %in% c("25", "50", "75"))


pred_vb <- expand_grid(
  han_vb_25_75,
  age = seq(0, 26, 0.5)
) %>% 
  mutate(
    # age_est = (log(1 - (tl / l_inf)) / -k ) + t0
    tl = l_inf * (1 - exp(-k *(age - t0))) 
  )



ggplot() + 
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm),
             alpha = 0.45, size = 3) +
  geom_line(data = pred_length_wide_df, 
            aes(x = pred_age, y = m_50), linewidth = 1, colour = "#39808f") + 
  geom_ribbon(data = pred_length_wide_df,
              aes(ymin = m_2.5,
                  ymax = m_97.5,
                  x = pred_age, y = m_50),
              fill = "#39808f", alpha = 0.25) +
  geom_line(data = pred_vb, aes(x = age, y = tl, 
                                linetype = percentile), 
            linewidth = 1) + 
  scale_linetype_manual(
    name = "Percentile", 
    values = c(2:5)
  ) +
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(0, 1000, 100)) + 
  theme_classic(base_size = 15) + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)") 

