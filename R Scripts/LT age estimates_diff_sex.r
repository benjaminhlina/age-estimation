here::i_am("R Scripts/LT age estimates.r")

# bring in packages -----
library(broom)
library(boot)
library(car)
library(data.table)
library(dplyr)
library(FSA)
library(FSAdata)
library(FSAmisc)
library(ggplot2)
library(here)
library(lubridate)
library(readr)
library(xlsx)
library(openxlsx)
library(tibble)
# bring in QC age estimtes with lengths 

lt_ae_qc <- read_csv(here("data", 
                          "qc_collection_2020.csv"))



lt_ae_bd <- read_csv(here("data", 
                          "bd_gn_2019_2020.csv"))


lt_ae_qc$cr <- "qc" 
lt_ae_bd$cr <- "bd" 

lt_ae <- rbind(lt_ae_qc, lt_ae_bd)




glimpse(lt_ae)

lt_ae <- lt_ae %>% 
  arrange(age_est)

lt_ae$year <- dmy(lt_ae$date) %>% 
  year()



lt_7 <- lt_ae %>% 
  filter(age_est == 7)

lt_7






# slim down the data ----

lt_slim <- lt_ae %>% 
  select(date, year, species:wt_g, sex, age_est, cr)
lt_slim$age_est <- as.integer(lt_slim$age_est)

lt_slim <- lt_slim %>% 
  arrange(tl_mm, age_est)


glimpse(lt_slim)

# look at summary of the data -----
tapply(lt_slim$age_est, lt_slim$year, summary)

?vbFuns



mf.start <- lt_slim %>% 
  filter(sex == "m") %>% 
  vbStarts(tl_mm ~ age_est,data = .)


ggplot(data = lt_slim, aes(x = age_est, y = tl_mm)) + 
  geom_point(alpha = 0.45, size = 3, aes(colour = sex)) +
  
  geom_smooth(data = lt_slim %>% 
                filter(sex == "m"),
    method = "nls", se = FALSE, fullrange = TRUE, 
    method.args = list(formula = y ~ Linf * (1 - exp(-K * (x - t0))),
                       start = list(Linf = mf.start$Linf, 
                                    K = mf.start$K, 
                                    t0 = mf.start$t0))) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0, 26, 2)) +
  # scale_y_continuous(breaks = seq(50, 900, 50), limits = c(50,)) +
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")



# lt_slim <- lt_slim %>%
#   filter(cr == "qc")



# look at von bertlanffy growth curve summary of exquation ----
(vb <- vbFuns(param = "Typical"))
predict2 <- function(x) predict(x, newdata = data.frame(age_est = a))
a <- seq(-1, 25, by = 0.5)

sex <- unique(lt_slim$sex)
sex

lts <- lt_slim %>% 
  filter(sex %in% c("f", "m"))

sx <- unique(lts$sex)
sx
out <- list()
model_output <- list()
model_ci <- list()
models <- list()
View(dats)

# creeat for loop to create qc and bd model -----
for (i in 1:length(sx)) {
  
  dats <- lt_slim %>% 
    filter(sex %in% sx[i])
  
  sexs <- unique(dats$sex)
  
  dats <- dats %>% 
    filter(age_est <= 15)
    p <- ggplot(data = dats, aes(x = age_est, y = tl_mm)) +
    geom_point(alpha = 0.45) +
    # geom_line(aes(y = fitted(f.fit))) + 
    theme_classic() + 
    labs(x = "Estimated Age (Yr)", 
         y = "Total Length (mm)")
  
  print(p)
  ?walfordPlot
  walfordPlot(tl_mm ~ age_est, data = dats)
  # determine L inf, K which is growth rate and t0 which is age at which lenght = 0 
  (f.starts <- vbStarts(tl_mm ~ age_est, data = dats))
  # run von b curve 
  f.fit <- nls(tl_mm ~ vb(age_est, Linf, K, t0),
               data = dats, start = f.starts)
  
  
  
  
  
  
  
  vb_output <- broom::tidy(f.fit)
  vb_confit <- Confint(f.fit) %>%
    as.data.frame() %>% 
    rownames_to_column("parameter") %>% 
    as_tibble()
  vb_output$sexs <- sexs
  vb_confit$sexs <- sexs
  
  
  
  
  # look at conffienceents, summary and CI for coefficencts 
  
  # use boot from Car to simulate nnumbers for model 
  
  
  # use function of predict within predict 2 functions 
  f.boot2 <- car::Boot(f.fit, f = predict2)
  
  preds1 <- data.frame(sexs,
                       a,
                       predict(f.fit, data.frame(age_est = a)),
                       Confint(f.boot2)) 
  
  
  names(preds1) <- c("sexs", "age", "fit", "LCI", "UCI")
  
  preds1 <- as_tibble(preds1)
  
  out[[i]] <- preds1
  model_output[[i]] <- vb_output
  model_ci[[i]] <- vb_confit
  models[[i]] <- f.fit
}


out_merge <- do.call(rbind, out)

out_merge
mo_merge <- do.call(rbind, model_output)
mo_merge

mci_mege <- do.call(rbind, model_ci)

mci_mege






# preds2 <- out_merge %>% 
#   filter(age >= 3 & age <= 25)

preds2 <- out_merge %>% 
  filter(case_when(crew %in% "bd" ~ age >= 5 & age <= 25 ,
                   crew %in% "qc" ~ age >= 3 & age <= 25))


#   geom_point(alpha = 0.45,) +
#   geom_line(aes(y = fitted(f.fit))) + 
#   theme_classic() + 
#   labs(x = "Estimated Age (Yr)", 
#        y = "Total Length (mm)")

p1 <- ggplot() +
  geom_ribbon(data = preds2, aes(x = age, ymin = LCI, ymax = UCI, 
                                 group = crew, fill = crew), alpha = 0.2) + 
  geom_point(data = lt_slim, aes(x = age_est, y = tl_mm, colour = cr), 
             alpha = 0.25, size = 3) + 
  geom_line(data = preds2, aes(x = age, y = fit, group = crew,
                               colour = crew), size = 1) + 
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



p1 + 
  annotate(geom = "text",
           label = makeVBEqnLabel(models[[1]]),
           parse = TRUE,
           size = 4,
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -1.5) + 
  annotate(geom = "text",
           label = makeVBEqnLabel(models[[2]]),
           parse = TRUE,
           size = 4,
           x = Inf, y = -Inf,
           hjust = 1.1, vjust = -0.5)


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



p <- ggplot(data = lts, aes(x = age_est, y = tl_mm, group = sex)) +
  geom_point(alpha = 0.45, aes(colour = sex), size = 4) +
  # geom_line(aes(y = fitted(f.fit))) + 
  theme_classic() + 
  labs(x = "Estimated Age (Yr)", 
       y = "Total Length (mm)")
p
