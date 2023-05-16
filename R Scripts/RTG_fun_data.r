library(dplyr)
library(ggplot2)
library(googlesheets4)
library(readr)
dat <- read_csv(here::here("RTG_Fish_list.csv")) %>%
  janitor::clean_names()
glimpse(dat)
dat <- dat %>%
  mutate(rtg = as.numeric(rtg))
ggplot(data = dat, aes(x = country_origin, y = rtg)) +
  geom_point(aes(colour = type_of_fish), size = 3) +
  scale_colour_viridis_d(end = 0.85) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90)