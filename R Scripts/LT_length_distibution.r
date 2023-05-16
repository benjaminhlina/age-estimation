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



write_rds(x = lt_slim, 
          file = here("Saved Data", 
                      "cleaned_length_at_age_raw.rds"))
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

lt_tot$basin <- factor(lt_tot$basin, levels = c("East Basin", 
                                                "West Basin", 
                                                "North Basin", 
                                                NA))
  
# View(lt_tot)

# plots ---------------------------------------------------------
ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), 
                 position = "identity", 
                 alpha = 0.5) + 
  facet_wrap(.~ year) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p 
p

ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), 
                 position = "identity", 
                 alpha = 0.5) + 
  # facet_wrap(.~ year) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p1
p1


ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), 
                 position = "identity", 
                 alpha = 0.5) + 
  facet_rep_wrap(.~ basin, repeat.tick.labels = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p2

p2

ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), colour = "black",
                 alpha = 0.5, position = "identity") + 

  facet_rep_wrap(.~ basin, repeat.tick.labels = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)
                     ) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p3

p3


ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), colour = "black", 
                 position = "stack", 
                 alpha = 0.5) + 
  facet_wrap(.~ year) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p4
p4

ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 35, aes(fill = year, group = year), colour = "black",
                 position = "stack", 
                 alpha = 0.5) + 
  # facet_wrap(.~ year) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 36, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 25)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 90)) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p5
p5


ggplot(data = lt_tot, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = year, group = year), colour = "black",
                 position = "stack", 
                 alpha = 0.5) + 
  facet_rep_wrap(.~ basin, repeat.tick.labels = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 26, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 25)) + 
  scale_fill_viridis_d(name = "Year", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 90)) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p6

p6

lt_tot_na <- lt_tot %>% 
  tidyr::drop_na()

ggplot(data = lt_tot_na, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, aes(fill = basin), 
                 colour = "black",
                 position = "stack", 
                 alpha = 0.75) + 
  facet_wrap(.~ year) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 36, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 25)) + 
  scale_fill_viridis_d(name = "Basin", option = "B", begin = 0.2, end = 0.7) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90)) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p7
p7



# ggsave(filename = here("Plots", "length_histogram_facet.png"), plot = p,
#        height = 10, width = 15)
# ggsave(filename = here("Plots", "length_histogram.png"), plot = p1,
#        height = 8.5, width = 11)
# ggsave(filename = here("Plots", "length_histogram_facet_basin.png"), plot = p2,
#        height = 10, width = 15)
# ggsave(filename = here("Plots", "tl_hist_facet_basin.png"), plot = p3,
#        height = 10, width = 15)
# ggsave(filename = here("Plots", "tl_facet_year_outline.png"), plot = p4,
#        height = 12, width = 15)
# ggsave(filename = here("Plots", "length_histogram_stack.png"), plot = p5,
#        height = 8.5, width = 11)
ggsave(filename = here("Plots", "tl_hist_facet_basin_stack.png"), plot = p6,
       height = 10, width = 15)
# ggsave(filename = here("Plots", "tl_facet_year_basin_outline.png"), plot = p7,
#        height = 12, width = 15)
