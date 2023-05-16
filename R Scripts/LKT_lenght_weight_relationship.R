here::i_am("R Scripts/LT age estimates.r")

# bring in packages -----
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
library(readr)
library(xlsx)
library(openxlsx)
library(tibble)
# bring in QC age estimtes with lengths 


fish_tag_data <- read_csv(here::here("Data", 
                                     "all fish tagged kenauk.csv")) %>%
  janitor::clean_names()

glimpse(fish_tag_data)

fish_tag_data <- fish_tag_data %>% 
  mutate(tag_date = dmy(tag_date), 
         year = year(tag_date))

#  grab just lake trout to add weights -----

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

# slim down the data ----

lt_slim <- lt_ae %>% 
  select(basin:wt_g, sex, age_est, cr)
lt_slim$age_est <- as.integer(lt_slim$age_est)

glimpse(lt_2021)





lt_slim <- lt_slim %>% 
  arrange(tl_mm, age_est)

lt_slim <- lt_slim %>% 
  select(tag_no, basin, tl_mm, fl_mm, girth_mm, wt_g) %>% 
  mutate(year = 2020)

glimpse(lt_2021)
lt_2021_s <- lt_2021 %>% 
  select(basin:comments)
lt_2021_s <- lt_2021_s %>% 
  select(tag_id, basin, tl_mm, fl_mm, girth_mm, wt_g) %>% 
  rename(tag_no = tag_id) %>% 
  mutate(year = 2021)


lt_tot <- bind_rows(lt_slim, lt_2021_s)
lt_tot$year <- as.factor(lt_tot$year)
lt_tot <- lt_tot %>% 
  arrange(tl_mm)

glimpse(lt_tot)

lt_tot <- lt_tot %>% 
  mutate(
    basin = factor(stringr::str_replace(as.character(basin), " Basin",
                                        ""), 
                   levels = c("East", "West", "North")),
    vemco_type = NA
  )

glimpse(lt_tot)
lt_tot <- lt_tot %>% 
  bind_rows(fish_acel)

levels(lt_tot$year)
glimpse(lt_tot)


# ----- log weight and lenght and regressioin ------

lt_tot <- lt_tot %>% 
  mutate(
    tl_log = log10(tl_mm), 
    wt_log = log10(wt_g), 
    year = forcats::fct_relevel(year, c("2017", "2018", "2019", 
                   "2020", "2021"))
  )


lt_tot_v <- lt_tot %>% 
  filter(year %in% c("2017", "2018", "2019"))


m <- lm(wt_log ~ tl_log, data = lt_tot)


car::Anova(m)
anova(m)
summary(m)

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


glimpse(predicts)
    
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
           "LKT_length_weight_regression.png"), 
      plot = p, 
      width = 11, 
      height = 8.5)


coef(m)
