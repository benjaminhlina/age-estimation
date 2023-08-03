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
  library(openxlsx)
  library(patchwork)
  library(purrr)
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
# m2 <- lm(wt_log ~ tl_log * year, data = lt_tot)


car::Anova(m)
anova(m)
summary(m)

mod.emt <- emtrends(m, ~3, var = "tl_log")
mod.emt
sig_3 <- test(mod.emt, null = 3, side = ">")

sig_3
sig_3 %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  write.xlsx(., here("Results", 
                     "tl_wt_allometric_vs_isometric_results.xlsx"))
# car::Anova(m2)
# anova(m2)
# summary(m2)
# 
# 
# m1 <- lmer(wt_log ~ tl_log + (1|year), data = lt_tot)
# 
# car::Anova(m1, type = "III")
# anova(m1)
# summary(m1)

# ---- model outputs -----

glance(m)
tidy(car::Anova(m))
tidy(m)

summary(m)
# 
# glance(m1)
# tidy(car::Anova(m))
# 
# tidy(m1)
# 
# summary(m1)


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
  tl_mm = seq(0, 770, 1)
) %>% 
  mutate(
    log_tl = log10(tl_mm), 
    log_w = log_a + (b * log_tl), 
    wt = 10 ^ log_w
  )



han_pred_max <- han_pred %>% 
  filter(wt < 5250 & 
           percentile %in% c("25", "50", "75"))

# ---- compare percentiles -----
pap_coef <- tibble(
  n = 1, 
  percentile = "Lake Papineau", 
  log_a = coef(m)[1],
  b = coef(m)[2]
)

comb_coef <- bind_rows(han_lw_wide, pap_coef)
comb_coef

# ---- From Hansen et al. 2021 we can calculate Kn ----
# to represent overall health 

# first grab the mean or 50 percentile 

han_pred_50 <- han_pred %>% 
  filter(percentile == "50")

han_pred

lt_tot_na <- lt_tot %>% 
  filter(wt_g != is.na(wt_g))

# join with lt_tot by tl_mm, han_pred has to be from every mm 
dat_k <- lt_tot_na %>% 
  left_join(han_pred_50, by = "tl_mm")

# calculate Kn by dividing our observed weights to 50% predicted weights 

dat_k <- dat_k %>% 
  mutate(
    kn = wt_g / wt
  ) %>% 
  filter(basin != is.na(basin))

# ---- evaulate kn disturbion and calculate general stats ---- 
ggplot(data = dat_k, aes(x = kn)) + 
  geom_histogram() + 
  facet_wrap(.~ percentile, scales = "free_x")

dat_k_50 <- dat_k %>% 
  filter(percentile %in% "50")
dat_k_50
# ---- create model to assess if Kn changes with length ----- 
fitdistrplus::descdist(dat_k_50$kn)

t <- fitdistrplus::fitdist(data = dat_k_50$kn, distr = "gamma",
                           method = "mme")
ts <- fitdistrplus::fitdist(data = dat_k_50$kn, distr = "norm", 
                            # method = "mme"
)

plot(t)
plot(ts)

# ---- model using GAMMA ---- 
m5 <- glm(kn ~ tl_mm * basin, data = dat_k_50,
          contrasts = list(basin = "contr.sum"),
          family = Gamma(link = "identity"), )


m6 <- update(m5, .~ tl_mm - basin)
m7 <- update(m5, .~ basin )

Anova(m7, type = "III")



lt_tot %>% 
  filter(wt_g != is.na(wt_g))
# create model list for model selection ------
model_list <- list(m5, m6, m7
)
# give the elements useful names
names(model_list) <- c("m5", "m6", "m7")


# get the summaries using `lapply

summary_list <- lapply(model_list, function(x) tidy(x, parametric = TRUE))
glance_list <- lapply(model_list, glance)

glance_summary <- map_df(glance_list, ~as.data.frame(.x), .id = "id") %>% 
  mutate(model = lapply(model_list, formula) %>%
           as.character() 
  ) %>% 
  dplyr::select(model, id:df.residual) %>% 
  arrange(AIC) %>% 
  mutate( 
         delta_AIC = AIC - first(AIC), 
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)


glance_summary


par(mfrow = c(2, 2))
plot(m6)
par(mfrow = c(1, 1))
res <- DHARMa::simulateResiduals(m6)
plot(residuals(m6))
hist(residuals(m6))
plot(res)


car::Anova(m6, type = "III")
summary(m6)




# ------------------ ONESIDED TEST OF IF THE RATIO IS GREATER THAN 1 -------


t.test(dat_k_50$kn, mu = 1, alternative = "less") 
kn_test <- tapply(dat_k_50$kn, dat_k_50$basin, function(x) t.test(x, mu = 1, alternative = "less")) 


kn_test_sum <- kn_test %>% 
  imap(~ broom::tidy(.x)) %>% 
  bind_rows(.id = "basin")
kn_test_sum

dat_k_50 %>% 
  group_by(
    basin
  ) %>% 
  summarise(
    kn_m = mean(kn),
    sem = sd(kn) / sqrt(n())
  ) %>% 
  ungroup() %>% 
  write.xlsx(., here("Results", 
                     "kn_mean.xlsx"))

write.xlsx(kn_test_sum, here("Results", 
                             "kn_t_test_basin.xlsx"))
# Kn changes with length 

# ---- summary table of Kn ---- 

kn_sum <- dat_k_50 %>% 
  filter(basin != is.na(basin) & kn != is.na(basin)) %>% 
  # group_by(
  #   basin
  # ) %>% 
  summarise(
    
    mean_kn = mean(kn), 
    sd_kn = sd(kn), 
    sem_kn = sd(kn) / sqrt(n()),
    n = n(), 
  )


kn_sum


dat_k_50 %>% 
  filter(kn > 1)

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



# ---- plot ------ 
predicts <- predicts %>% 
  filter(basin != is.na(basin))

ggplot(data = predicts) + 
  geom_point(aes(x = tl_log, y = wt_log, 
                 colour = basin), 
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
  ) -> p1

p1

ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_regression_basin.png"), 
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

glimpse(han_pred_max)

han_pred_max <- han_pred_max %>% 
  mutate(
    percentile = forcats::fct_rev(percentile)
  )

# ---- plot back transformed ----- 
ggplot(data = predicts) + 
  geom_point(aes(x = tl_mm, y = wt_g, 
                 colour = basin), 
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
    linewidth = 1,
    aes(x = tl_mm, y = wt, linetype = percentile)
  ) + 
  scale_colour_viridis_d(name = "Basin", 
                         begin = 0.25, end = 0.75, option = "D", 
                         # alpha = 0.7
                         alpha = 0.7
  ) + 
  scale_y_continuous(breaks = seq(0, 5000, 500)) +
  scale_x_continuous(breaks = seq(100, 700, 100)) + 
  theme_bw(base_size = 15) + 
  theme(
    axis.text = element_text(colour = "black"), 
    legend.position = c(0.07, 0.65), 
    legend.background = element_blank(), 
    panel.grid = element_blank(),
    plot.tag.position  = c(0.11, 0.97)
  ) + 
  scale_linetype_manual(
    name = "Percentile", 
    values = c(2:5)
  ) +
  ggpmisc::stat_poly_eq(aes(x = tl_log, y = wt_log,
                            label = paste(after_stat(eq.label), 
                                          after_stat(rr.label), sep = "*\", \"*")),
                        label.x = 0.95,
                        label.y = 0.03) +
  labs(
    x = "Total Length (mm)", 
    y = "Weight  (g)"
  ) -> p2
p2


ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_non_transform_regress.png"), 
       plot = p2, 
       width = 11, 
       height = 8.5)


# ---- plot Kn vs length ----- 
1 - (m6$deviance/m6$null.deviance)
# ----- create - predicted dataframe ------


ggplot(
  data = dat_k_50, aes(y = kn, x = tl_mm)
) +
  geom_hline(yintercept = 1, linetype = 2) + 
  stat_smooth(method = "lm", colour = "black", se = TRUE) +
  geom_point(data = dat_k_50, aes(y = kn, x = tl_mm, colour = basin), 
             size = 3
  ) + 
  scale_colour_viridis_d(name = "Basin", 
                         begin = 0.25, end = 0.75, alpha = 0.7) +
  scale_x_continuous(breaks = seq(100, 700, 100)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.tag.position  = c(0.11, 0.97)
  ) +
  ggpmisc::stat_poly_eq(aes(x = tl_mm, y = kn,
                            label = paste(after_stat(eq.label),
                                          after_stat(rr.label), sep = "*\", \"*")),
                        label.x = 0.98,
                        label.y = 0.03) +
  labs(
    x = "Total Length (mm)", 
    y = expression(K[n]) 
  ) -> p3 

p3
ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_Kn.png"), 
       plot = p3, 
       width = 11, 
       height = 8.5)




ggplot(data = dat_k, aes(y = kn, x = tl_mm)) +
  # geom_hline(yintercept = 1, linetype = 2) + 
  # stat_smooth(method = "lm") + 
  geom_point(size = 4, aes(colour = basin)
  ) + 
  scale_colour_viridis_d(name = "Basin", 
                         begin = 0.35, end = 0.75, alpha = 0.7) +
  scale_x_continuous(breaks = seq(100, 700, 100)) + 
  facet_wrap(.~ percentile, scales = "free_y") + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        legend.position = "none"
  ) +
  labs(
    x = "Total Length (mm)", 
    y = expression(K[n]) 
  ) -> p4

# p4
ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_Kn_percentile.png"), 
       plot = p4, 
       width = 11 * 1.25, 
       height = 8.5)



unique(han$metric)


ggplot(data = dat_k_50, aes(x = basin, y = kn)) + 
  geom_jitter(aes(colour = basin), width = 0.1, size = 3) +
  geom_hline(yintercept = 1, linetype = 2) + 
  scale_colour_viridis_d(name = "Basin", 
                         begin = 0.25, end = 0.75, alpha = 0.7) +
  # scale_x_continuous(breaks = seq(100, 700, 100)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.tag.position  = c(0.11, 0.97)
  ) +
  labs(
    x = "Capture Basin", 
    y = expression(K[n]) 
  ) -> p6

p6


p5 <- p2 / p3 + 
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")
p7 <- p2 / p6 + 
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")

# p5
ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_model_w_kn.png"), 
       width = 8.5, height = 11, units = "in", plot = p5, 
       dpi = 300)
ggsave(here("Plots",
            "length weight relationship",
            "LKT_length_weight_model_w_kn_group.png"), 
       width = 8.5, height = 11, units = "in", plot = p7, 
       dpi = 300)


