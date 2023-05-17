# ---- bring in packages -----
{
  library(broom)
  library(boot)
  library(car)
  library(data.table)
  library(dplyr)
  library(FSA)
  library(FSAdata)
  library(ggplot2)
  library(here)
  library(lubridate)
  library(lemon)
  library(openxlsx)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
}

# ---- bring data ---- 

# bring table from Hanset et al. 2021 in Muir et al. 2021. 

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

# bring in QC age estimates with lengths 


fish_tag_data <- read_csv(here::here("Data", 
                                     "all fish tagged kenauk.csv")) %>%
  janitor::clean_names()

glimpse(fish_tag_data)

fish_tag_data <- fish_tag_data %>% 
  mutate(tag_date = dmy(tag_date), 
         year = year(tag_date))

# ---- grab just lake trout to add weights -----

fish_acel <- fish_tag_data %>%
  filter(species %in% "LT") %>% 
  rename(tag_no = gray_floy_tag_number) %>% 
  dplyr::select(tag_no, basin, tl, fl, girth, weight, year, vemco_type) %>% 
  rename(
    tl_mm = tl, 
    fl_mm = fl,
    girth_mm = girth,
    wt_g = weight
  ) %>% 
  mutate(
    basin = factor(stringr::str_replace_all(as.character(basin), 
                                            c(" Basin" = "", "Main" = "East")),
                   levels = c("East", "West", "North")), 
    year = as.factor(year)
    
  )


glimpse(fish_acel)

# ---- bring in gn sampled fish ----- 
lt_ae_qc <- read_csv(here("Data", 
                          "qc_collection_2020.csv"))



lt_ae_bd <- read_csv(here("Data", 
                          "bd_gn_2019_2020.csv"))


lt_2021 <- read_csv(here("Data", 
                         "gill_net_fish_data_2021.csv"))

lt_ae_qc$cr <- "qc" 
lt_ae_bd$cr <- "bd" 


glimpse(lt_ae_qc)
glimpse(lt_ae_bd)

lt_ae <- rbind(lt_ae_qc, lt_ae_bd)

glimpse(lt_ae)

# ---- slim down the data ----

lt_slim <- lt_ae %>% 
  select(basin:wt_g, sex, age_est, cr)
# make age as integer 
lt_slim$age_est <- as.integer(lt_slim$age_est)


lt_slim <- lt_slim %>% 
  arrange(tl_mm, age_est)

# arragne so tag number is in front of basin add in year 
lt_slim <- lt_slim %>% 
  select(tag_no, basin, tl_mm, fl_mm, girth_mm, wt_g) %>% 
  mutate(year = 2020)

glimpse(lt_2021)
# select columns we need from 2021 to join with 2020 collection 
lt_2021_s <- lt_2021 %>% 
  select(basin:comments)

lt_2021_s <- lt_2021_s %>% 
  select(tag_id, basin, tl_mm, fl_mm, girth_mm, wt_g) %>% 
  rename(tag_no = tag_id) %>% 
  mutate(year = 2021)


# ---- bind gn lt 2020 and 2021
lt_tot <- bind_rows(lt_slim, lt_2021_s)
lt_tot$year <- as.factor(lt_tot$year)

lt_tot <- lt_tot %>% 
  arrange(tl_mm)

glimpse(lt_tot)

# fix basin namining 
lt_tot <- lt_tot %>% 
  mutate(
    basin = factor(stringr::str_replace(as.character(basin), " Basin",
                                        ""), 
                   levels = c("East", "West", "North")),
    vemco_type = NA
  )

glimpse(lt_tot)
# ---- binding 2017-2019 data with 2020 and 2021 data ---- 
lt_tot <- lt_tot %>% 
  bind_rows(fish_acel)

levels(lt_tot$year)
glimpse(lt_tot)


# ---- log weight and length and regression ------

lt_tot <- lt_tot %>% 
  mutate(
    tl_log = log10(tl_mm), 
    wt_log = log10(wt_g), 
    year = forcats::fct_relevel(year, c("2017", "2018", "2019", 
                                        "2020", "2021"))
  )


lt_tot_v <- lt_tot %>% 
  filter(year %in% c("2017", "2018", "2019"))

# ---- run the regression ----- 

m <- lm(wt_log ~ tl_log, data = lt_tot)


car::Anova(m)
anova(m)
summary(m)

# ---- model outputs -----

glance(m)
tidy(car::Anova(m))





# ---- predict -----


dat_2 <- lt_tot

fits <- predict(m, newdata = dat_2, 
                type = "response", se.fit = TRUE)

# combine fits with dataframe for plotting and calc upper and lower 
# add in month abb for plotting 
predicts <- data.frame(dat_2, fits) %>% 
  mutate(
    lower = fit - 1.96 * se.fit,
    upper = fit + 1.96 * se.fit
  )

# ---- use Hansen et al. parameters to estimate quantiles ---- 
glimpse(han)

unique(han$metric)

# grab only metrics we need 

han_lw <- han %>% 
  filter(metric %in% c("log10(α)", "β"))

han_lw_wide <- han_lw %>% 
  pivot_wider(
    names_from = metric, values_from = value
  ) %>% 
  rename(
    log_a = `log10(α)`,
    b = β
  ) %>% 
  mutate(
    log_a = (as.numeric(gsub("\\–", "", log_a)) * -1),
    b = as.numeric(b)
    )

han_lw_wide
han_pred <- expand_grid(
  han_lw_wide,
  tl = seq(0, 770, 10)
) %>% 
  mutate(
    log_tl = log10(tl), 
    log_w = log_a + (b * log_tl), 
    wt = 10 ^ log_w
  )

han_pred_max <- han_pred %>% 
  filter(wt < 5250 & 
           percentile %in% c("25", "50", "75", "97.5"))
# ---- plot ------ 
ggplot(data = predicts) + 
  geom_point(aes(x = tl_log, y = wt_log, 
                 colour = year), 
             size = 3) + 
  # stat_smooth(geom = "line", method = "lm",
  #             aes(x = tl_log, y = wt_log),
  #             colour = "black",
  #             linewidth = 1, se = TRUE) +
  geom_line(colour = "black", linewidth = 1, aes(x = tl_log, y = fit)) +
  geom_ribbon( 
    aes(ymin = lower,
        ymax = upper,
        x = tl_log, y = fit), alpha = 0.25) +
  scale_colour_viridis_d(name = "Year", 
                         begin = 0.17, end = 0.8, alpha = 0.7) + 
  theme_classic(base_size = 15) + 
  theme(
    axis.text = element_text(colour = "black"), 
    legend.position = c(0.06, 0.90), 
    legend.background = element_blank()
  ) + 
  ggpmisc::stat_poly_eq(aes(x = tl_log, y = wt_log,
                            label = paste(after_stat(eq.label),
                                          after_stat(rr.label), sep = "*\", \"*")), 
                        label.x = 0.03, 
                        label.y = 0.80) +
  labs(
    x = expression(paste(Log[10], "[Total Length (mm)]" )), 
    y = expression(paste(Log[10], "[Weight  (g)]" ))
  ) -> p 

p

ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_regression.png"), 
       plot = p, 
       width = 11, 
       height = 8.5)




glimpse(predicts)
# ---- back transform fit and confidence intervals 
predicts <- predicts %>% 
  mutate(
    anti_fit = 10 ^ fit,
    anti_lower = 10 ^ lower,
    anti_upper = 10 ^ upper
    
  )
# predict_2 <- predict_2 %>% 
#   mutate(
#     anti_fit = 10 ^ fit,
#     anti_lower = 10 ^ lwr,
#     anti_upper = 10 ^ upr
#     
#   )

glimpse(lt_tot)

# ---- plot back transformed ----- 
ggplot(data = predicts) + 
  geom_point(aes(x = tl_mm, y = wt_g, 
                 colour = year), 
             size = 3) + 
  # stat_smooth(geom = "line", method = "lm",
  #             aes(x = tl_log, y = wt_log),
  #             colour = "black",
  #             linewidth = 1, se = TRUE) +
  geom_line(colour = "black", linewidth = 1, aes(x = tl_mm, y = anti_fit)) +
  geom_ribbon(
    aes(ymin = anti_lower,
        ymax = anti_upper,
        x = tl_mm, y = anti_fit), alpha = 0.25) +
  geom_line(
    data = han_pred_max, 
    aes(x = tl, y = wt, linetype = percentile)
  ) + 
  scale_colour_viridis_d(name = "Year", 
                         begin = 0.17, end = 0.8, alpha = 0.7) + 
  scale_y_continuous(breaks = seq(0, 5000, 500)) +
  scale_x_continuous(breaks = seq(100, 700, 100)) + 
  theme_classic(base_size = 15) + 
  theme(
    axis.text = element_text(colour = "black"), 
    legend.position = c(0.06, 0.8), 
    legend.background = element_blank()
  ) + 
  scale_linetype_manual(
    name = "Percentile", 
    values = c(2:5)
  ) +
  ggpmisc::stat_poly_eq(aes(x = tl_log, y = wt_log,
                            label = paste(after_stat(eq.label),
                                          after_stat(rr.label), sep = "*\", \"*")),
                        label.x = 0.02,
                        label.y = 0.80) +
  labs(
    x = "Total Length (mm)", 
    y = "Weight  (g)"
  ) -> p2
# p2


ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_non_transform_regress.png"), 
       plot = p2, 
       width = 11, 
       height = 8.5)
