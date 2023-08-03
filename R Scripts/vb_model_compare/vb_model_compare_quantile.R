# ---- bring in packages -----
{
  library(broom)
  library(broom.mixed)
  library(car)
  library(data.table)
  library(dplyr)
  library(emmeans)
  library(FSA)
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

lt_slim <- read_rds(here("Saved Data", 
                         "cleaned_length_at_age_raw.rds"))

p2 <- read_rds(here("Saved Plots", 
                    "posterior distribution_lkt_age_length.rds"))

glimpse(lt_slim)

lkt_s <- lt_slim %>% 
  select(age_est, tl_mm) %>% 
  rename(age = age_est, 
         tl = tl_mm) %>% 
  mutate(
    lake_name = "Papineau", 
    linf = NA, 
    t0 = NA, 
    k = NA
  ) %>% 
  select(lake_name, linf:k, age, tl)
# ---- select vb curve parameters and remove rows with NA -----
vb <- han_raw %>% 
  select(lake_name, linf, t0, k)

vb_han <- vb %>% 
  drop_na()


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
vb_all <- bind_rows(vb_han, 
                    vb_pl_coef)

vb_all_p <- expand_grid(
  vb_all,
  age = seq(0, 26, 0.1)
) %>% 
  mutate(
    tl = linf * (1 - exp(-k *(age - t0))) 
  )

pred_vb <- vb_all_p %>% 
  filter(lake_name != "Papineau")

# ----- plot preliminary VB to see where pap comes in ----- 

ggplot() + 
  geom_point(data = pred_vb, aes(x = age, y = tl, 
                                 # group = lake_name,
                                 # colour = lake_name
  ), size = 3) + 
  geom_line(data = vb_pl_wide, aes(x = pred_age, y = m_50), 
            colour = "black", linewidth = 1, linetype = 2) + 
  scale_colour_viridis_d(end = 0.85, name = "Lakes") + 
  geom_quantile(data = pred_vb, aes(x = age, y = tl), colour = "red",
                quantiles = seq(0.25, 0.75, 0.25), method = "rqss", 
                lambda = 1, linewidth = 1)




# ----- quantile regression -----
(vb <- vbFuns(param = "Typical"))

(f.starts <- vbStarts(tl ~ age, data = vb_all_p))

qr <- c(0.25, 0.50, 0.75)
qr_list <- list()

for (i in 1:length(qr)) {
  
  vb_qr <- nlrq(tl ~ vb(age, linf, k, t0),
                tau = 0.25, start = list(
                  linf = 776.186, 
                  k = 0.0124, 
                  t0 = -8.846
                ), 
                data = vb_all_p, trace = TRUE)
  
  qr_list[[i]] <- vb_qr
}
anova(qr_list, test = "Wald", joint=TRUE)
car::Anova(qr_list[[1]], qr_list[[2]])
AIC(qr_list[[1]])
AIC(qr_list[[2]])
AIC(qr_list[[3]])
predict
predict(qr_list[[1]])

qr_list %>% 
  imap(~ broom.mixed::glance(.x))
anova.(vb_qr)
qr_est <- qr_list %>% 
  imap(~ broom.mixed::tidy(.x)) %>% 
  bind_rows(.id = "quantile") %>% 
  mutate(
    quantile = factor(case_when(
      quantile == "1" ~ "25",
      quantile == "2" ~ "50",
      quantile == "3" ~ "75",
    ), levels = c("75", "50", "25"))
  )

compare_qr <- vb_pl_coef %>%
  pivot_longer(
    -lake_name, 
    names_to = "term", 
    values_to = "estimate"
  ) %>% 
  rename(
    quantile = lake_name
  ) %>% 
  mutate(
    std.error = NA, 
    statistic = NA,
    p.value = NA
  )



com_qr <- bind_rows(qr_est, 
                    compare_qr)


qr_pred <- qr_est %>% 
  pivot_wider(
    id_cols =  quantile,
    names_from = "term", 
    values_from = c(estimate)) %>% 
  expand_grid(
    age = seq(0, 26, 0.1)
  ) %>% 
  mutate(
    tl = linf * (1 - exp(-k *(age - t0))) 
  )

qr_predict <- qr_list %>% 
  imap(~broom.mixed::augment(.x)) %>% 
  bind_rows(.id = "quantile") %>% 
  mutate(
    quantile = factor(case_when(
      quantile == "1" ~ "25",
      quantile == "2" ~ "50",
      quantile == "3" ~ "75",
    ), levels = c("75", "50", "25"))
  )


qr_pred_pap <- vb_pl_wide %>%
  select(pred_age, m_50) %>% 
  rename(age = pred_age, 
         tl = m_50) %>% 
  mutate(
    quantile = "Papineau", 
    linf = NA, 
    k = NA, 
    t0 = NA
  ) %>% 
  select(quantile:t0, age, tl)

qr_pred_pap <- bind_rows(qr_pred, 
                         qr_pred_pap)


ggplot(data = qr_pred_pap, aes(x = tl)) + 
  geom_histogram()

fitdistrplus::descdist(qr_pred_pap$tl)

library(glmmTMB)
library(DHARMa)

m <- glmmTMB(tl ~ age * quantile, data = qr_pred_pap, 
             family = t_family())


hist(resid(m))

res <- simulateResiduals(m)
plot(res)


multi_comp <- emmeans(m, ~ age * quantile, 
                      adjust = "Tukey")
contrast_effects <- contrast(multi_comp, method = "pairwise", adjust = "Tukey")
tidy(contrast_effects) %>% 
  janitor::clean_names() %>% 
  arrange(adj_p_value, contrast)  %>% 
  
  write.xlsx(., here("results", 
                     "vb_quantile_pap_compare_han_constrast.xlsx"))


# extract each parmater indivually to put von B. curve exquation on plot
linf <- formatC(vb_pl_coef$linf, format="f", digits = 0)
k <- formatC(vb_pl_coef$k, format = "f", digits = 3)
# Handle t0 differently because of minus in the equation
t0 <- vb_pl_coef$t0
t0 <- paste0(ifelse(t0 < 0, "+", "-"), formatC(abs(t0), format = "f", 
                                               digits = 3))
# Put together and return
labels <- paste("TL ==", linf,"~bgroup('(',1-e^{-", k,"~(age", t0,")},')')")

ggplot() +
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm),
             alpha = 0.45, size = 3) +
  geom_line(data = qr_pred, aes(x = age, y = tl, linetype = quantile), 
            linewidth = 1) + 
  # geom_line(data = qr_predict, aes(x = age, y = .fitted, linetype = quantile)) +
  geom_line(data = vb_pl_wide, 
            aes(x = pred_age, y = m_50), linewidth = 1, colour = "#39808f") + 
  geom_ribbon(data = vb_pl_wide,
              aes(ymin = m_2.5,
                  ymax = m_97.5,
                  x = pred_age, y = m_50),
              fill = "#39808f", alpha = 0.25) + 
  scale_linetype_manual(
    name = "Percentile", 
    values = c(2:5)
  ) +
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(0, 1000, 100)) + 
  annotate(geom = "text",
           label = labels,
           parse = TRUE,
           size = 4,
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -0.5) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = c(0.0725, 0.82), 
    panel.grid = element_blank(),
    plot.tag.position = c(0.12, 0.97),
    legend.background = element_blank()
    
  ) +
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)") -> p



p

# ----- plot update quantile   
dev.off()

png(filename =  here("Plots",
                     "von B Bayes",
                     "lt_bayes_vb_qr_with_posterior_stack.png"), 
    width = 8.5, height = 11, units = "in", 
    res = 300)

p6 <- p / p2 + 
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")

# p5
print(p6)
dev.off()



