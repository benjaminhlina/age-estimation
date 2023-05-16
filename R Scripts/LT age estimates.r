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
library(readr)
library(xlsx)
library(openxlsx)
library(tibble)
# bring in QC age estimtes with lengths 

lt_ae_qc <- read_csv(here("Data", 
                       "qc_collection_2020.csv"))



lt_ae_bd <- read_csv(here("Data", 
                          "bd_gn_2019_2020.csv"))


glimpse(lt_ae_qc)
lt_ae_qc$cr <- "qc" 
lt_ae_bd$cr <- "bd" 

lt_ae <- rbind(lt_ae_qc, lt_ae_bd)


tail(lt_ae_qc)

glimpse(lt_ae)

lt_ae <- lt_ae %>% 
  arrange(age_est)



lt_7 <- lt_ae %>% 
  filter(age_est == 7)

lt_7


# slim down the data ----

glimpse(lt_ae)

lt_ae %>% 
  ggplot(aes(x = tl_mm, fill = basin), colour = "black") +
  geom_histogram()

lt_slim <- lt_ae %>% 
  select(species:wt_g, basin, sex, age_est, cr)
lt_slim$age_est <- as.integer(lt_slim$age_est)

lt_slim <- lt_slim %>% 
  arrange(tl_mm, age_est)

lt_slim %>% 
  filter(age_est > 15)

glimpse(lt_slim)




write_rds(lt_slim, here("Saved Data", 
                        "cleaned_length_at_age_raw.rds"))

is.na(unique(lt_slim$wt_g))
summary(lt_slim)

lt_wt <- lt_slim %>% 
  filter(wt_g != is.na(wt_g))



ggplot(data = lt_slim, aes(x = tl_mm)) + 
  geom_histogram(bins = 30, fill = "black") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 14, 2)) + 
  scale_x_continuous(breaks = seq(150, 800, 50)) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Total Length (mm)", 
       y = "Frequency") -> p 

p



# look at summary of the data -----
tapply(lt_slim$age_est, lt_slim$sex, summary)

?vbFuns

# lt_slim <- lt_slim %>%
#   filter(cr == "qc")



# look at von bertlanffy growth curve summary of exquation ----
(vb <- vbFuns(param = "Typical"))

# determine L inf, K which is growth rate and t0 which is age at which lenght = 0 
(f.starts <- vbStarts(tl_mm ~ age_est, data = lt_slim))

?nls

# run von b curve 
f.fit <- nls(tl_mm ~ vb(age_est, Linf, K, t0),
             data = lt_slim, start = f.starts)


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) +
  geom_point(alpha = 0.45) +
  geom_line(aes(y = fitted(f.fit))) + 
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")




vb

vb_output <- broom::tidy(f.fit)
vb_confit <- Confint(f.fit) %>%
  as.data.frame() %>% 
  rownames_to_column("parameter") %>% 
  as_tibble()
vb_confit 
vb_output


# write.xlsx(vb_output, here("Results",
#                            "vb_model_output.xlsx"))
# write.xlsx(vb_confit, here("Results",
#                            "vb_ci_output.xlsx"))


# look at conffienceents, summary and CI for coefficencts 

# use boot from Car to simulate nnumbers for model 
f.boot1 <- Boot(f.fit)  # Be patient! Be aware of some non-convergence

vb_boot_coef <- Confint(f.boot1)
vb_boot_coef <- vb_boot_coef %>% 
  as_tibble() %>% 
  mutate_all(function(x) as.numeric(x)) %>% 
  mutate(metric = row.names(vb_boot_coef)) %>% 
  janitor::clean_names() %>% 
  rename(percent_2.5 = x2_5_percent, 
         percent_97.5 = x97_5_percent) %>% 
  select(metric, estimate:percent_97.5)

write.xlsx(vb_boot_coef, here("Results",
                           "vb_boot_ci_output.xlsx"))


# example of predcited values from models 
predict(f.fit, data.frame(age_est = 2:7))

# create function tobe used 
predict2 <- function(x) predict(x, newdata = data.frame(age_est = a))
a <- seq(-1, 25, by = 0.5)

# use function of predict within predict 2 functions 
f.boot2 <- car::Boot(f.fit, f = predict2)

preds1 <- data.frame(a,
                     predict(f.fit, data.frame(age_est = a)),
                     Confint(f.boot2)) 


names(preds1) <- c("age", "fit", "LCI", "UCI")

preds1 <- as_tibble(preds1)

preds1


preds2 <- preds1 %>% 
  filter(age >= 3 & age <= 25)



p1 <- ggplot() +
  geom_ribbon(data = preds2, aes(x = age, ymin = LCI, ymax = UCI),
              fill = "gray90") + 
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm), 
             alpha = 0.25, size = 3) + 
  geom_line(data = preds2, aes(x = age, y = fit), size = 1) + 
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(50, 900, 50)) + 
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")

p
p1

makeVBEqnLabel <- function(fit) {
  # Isolate coefficients (and control decimals)
  cfs <- coef(fit)
  Linf <- formatC(cfs[["Linf"]],format="f",digits=1)
  K <- formatC(cfs[["K"]],format="f",digits=3)
  # Handle t0 differently because of minus in the equation
  t0 <- cfs[["t0"]]
  t0 <- paste0(ifelse(t0<0,"+","-"),formatC(abs(t0),format="f",digits=3))
  # Put together and return
  paste0("TL ==",Linf,"~bgroup('(',1-e^{-",K,"~(age",t0,")},')')")
}



p2 <- p1 + 
  annotate(geom = "text",
           label = makeVBEqnLabel(f.fit),
           parse = TRUE,
                     size = 4,
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -0.5)



ggsave(here("Plots",
            "growth_curve_lt_papineau_t.png"), p2, width = 11, height = 8.5)
# ggsave(here("Plots", 
#             "growth_curve_lt_papineau.pdf"), p2, width = 11, height = 8.5)
ggplot(data = lt_slim, aes(x = tl_mm, y = wt_g)) + 
  geom_point() + 
  geom_smooth()





lt_slim <- lt_slim %>% 
  mutate(log_tl = log(tl_mm), 
         log_wt = log(wt_g))


m <- lm(data = lt_slim, log_tl ~ log_wt)

summary(m)


ggplot(data = lt_slim, aes(x = log_tl, y = log_wt)) + 
  geom_point(size = 3, colour = "blue") + 
  geom_smooth(method = "lm", colour = "black")


ggplot(data = lt_slim, aes(x = age_est)) + 
  geom_histogram(bins = 23, 
                 ) + 
  scale_x_continuous(breaks = seq(3, 23, 1))
hist(lt_slim$age_est)


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) + 
  geom_point(alpha = 0.45, size = 3) +
  
  geom_smooth( 
    method = "nls", se = FALSE, fullrange = TRUE, 
    method.args = list(formula = y ~ Linf * (1 - exp(-K * (x - t0))),
                       start = list(Linf = f.starts$Linf, 
                                    K = f.starts$K, 
                                    t0 = f.starts$t0))) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0, 26, 2)) +
  # scale_y_continuous(breaks = seq(50, 900, 50), limits = c(50,)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")

p3 <- ggplot() +
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm), 
             alpha = 0.25, size = 4) +  
  scale_x_continuous(breaks = seq(0, 26, 2)) + 
  scale_y_continuous(breaks = seq(50, 900, 50)) + 
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")
p3
