#library(devtools)
#devtools::install_github("wfinsinger/tapas")
library(mgcv)
library(tidyverse)
library(tapas) #https://github.com/wfinsinger/tapas/
source("mgcv2tapas_krzypl.R")

############tl27##############
tl27_sed <- read_csv("data/tl27_sed.csv")
tl27_age_tapas <- read_csv("data/tl27_age_tapas.csv")

tl27_tapas <- tl27_sed %>% 
  select(cmTop, cmBot, Sand) %>% 
  mutate(AgeTop = tl27_age_tapas$AgeTop, .after = cmBot) %>% 
  mutate(AgeBot = tl27_age_tapas$AgeBot, .after = AgeTop) %>% 
  mutate(Volume = 1, .after = AgeBot)

tl27_intercept <- median(diff(tl27_tapas$AgeTop))
tl27_pretreat <- pretreatment_data(tl27_tapas, yrInterp = tl27_intercept)

tl27_gam_pre <- tapas2mgcv(tl27_pretreat)
tl27_gam_pre <- transform(tl27_gam_pre, age = age_top - ((age_top-age_bot)/2))
tl27_gam_pre <- tl27_gam_pre[-length(tl27_gam_pre$cm_top), ] #remove the last row containing sample of thickness = 0
tl27_gam_pre$age <- round(tl27_gam_pre$age, digits = 3)
tl27_age_span <- max(tl27_gam_pre$age) - min(tl27_gam_pre$age)

set.seed(12)
tl27_gam <- gam(SandAR ~ s(age, k = 50), data = tl27_gam_pre, method = "REML")
tl27_gam_post <- mgcv2tapas_krzypl(tl27_gam) #mgcv2tapas function was modified as the original one gave en error (see file "mgcv2tapas_krzypl.R")  
summary(tl27_gam)

tl27_thresh_gl <- global_thresh(
  series = tl27_gam_post,
  proxy = "SandAR",
  thresh.value = 0.999,
  smoothing.yr = tl27_age_span,
  keep_consecutive = FALSE,
  minCountP = 0.05,
  MinCountP_window = tl27_intercept*5
)

Plot.Anomalies(
  series = tl27_thresh_gl,
  plot.crosses = TRUE,
  plot.x = TRUE,
  plot.neg = FALSE)

#########tl18-2###########
tl18_sed <- read_csv("data/tl18_sed.csv")
tl18_age_tapas <- read_csv("data/tl18_age_tapas.csv")

tl18_tapas <- tl18_sed %>% 
  select(cmTop, cmBot, Sand) %>% 
  mutate(AgeTop = tl18_age_tapas$AgeTop, .after = cmBot) %>% 
  mutate(AgeBot = tl18_age_tapas$AgeBot, .after = AgeTop) %>% 
  mutate(Volume = 1, .after = AgeBot)

tl18_intercept <- median(diff(tl18_tapas$AgeTop))
tl18_pretreat <- pretreatment_data(tl18_tapas, yrInterp = tl18_intercept)

tl18_gam_pre <- tapas2mgcv(tl18_pretreat)
tl18_gam_pre <- transform(tl18_gam_pre, age = age_top - ((age_top-age_bot)/2))
tl18_gam_pre <- tl18_gam_pre[-c(101:103), ] #delete the three last rows containing zeros in sandAR
tl18_gam_pre$age <- round(tl18_gam_pre$age, digits = 3)
tl18_age_span <- max(tl18_gam_pre$age) - min(tl18_gam_pre$age)

set.seed(12)
tl18_gam <- gam(SandAR ~ s(age, k = 100), data = tl18_gam_pre, method = "REML")
tl18_gam_post <- mgcv2tapas_krzypl(tl18_gam)  
summary(tl18_gam)

tl18_thresh_gl <- global_thresh(
  series = tl18_gam_post,
  proxy = "SandAR",
  thresh.value = 0.999,
  smoothing.yr = tl18_age_span,
  keep_consecutive = FALSE,
  minCountP = 0.05,
  MinCountP_window = tl18_intercept*5
)

Plot.Anomalies(
  series = tl18_thresh_gl,
  plot.crosses = TRUE,
  plot.x = TRUE,
  plot.neg = FALSE)

#####################tl09############
tl09_sed <- read_csv("data/tl09_sed.csv")
tl09_age_tapas <- read_csv("data/tl09_age_tapas.csv")

tl09_tapas <- tl09_sed %>% 
  select(cmTop, cmBot, Sand) %>% 
  mutate(AgeTop = tl09_age_tapas$AgeTop, .after = cmBot) %>% 
  mutate(AgeBot = tl09_age_tapas$AgeBot, .after = AgeTop) %>% 
  mutate(Volume = 1, .after = AgeBot)

tl09_tapas <- tl09_tapas[-c(32,33),] #remove observations corresponding to the bottom sand layer

tl09_intercept <- median(diff(tl09_tapas$AgeTop))

tl09_pretreat <- pretreatment_data(tl09_tapas, yrInterp = tl09_intercept)

tl09_gam_pre <- tapas2mgcv(tl09_pretreat)
tl09_gam_pre <- transform(tl09_gam_pre, age = age_top - ((age_top-age_bot)/2))
tl09_gam_pre$age <- round(tl09_gam_pre$age, digits = 3)
tl09_age_span <- max(tl09_gam_pre$age) - min(tl09_gam_pre$age)

set.seed(12)
tl09_gam <- gam(SandAR ~ s(age, k = 30), data = tl09_gam_pre, method = "REML")
tl09_gam_post <- mgcv2tapas_krzypl(tl09_gam)  
summary(tl09_gam)

tl09_thresh_gl <- global_thresh(
  series = tl09_gam_post,
  proxy = "SandAR",
  thresh.value = 0.999,
  smoothing.yr = tl09_age_span,
  keep_consecutive = FALSE,
  minCountP = 0.05,
  MinCountP_window = tl09_intercept*5
)

Plot.Anomalies(
  series = tl09_thresh_gl,
  plot.crosses = TRUE,
  plot.x = TRUE,
  plot.neg = FALSE)

#############tl08##############
tl08_sed <- read_csv("data/tl08_sed.csv")
tl08_age_tapas <- read_csv("data/tl08_age_tapas.csv")

tl08_tapas <- tl08_sed %>% 
  select(cmTop, cmBot, Sand) %>% 
  mutate(AgeTop = tl08_age_tapas$AgeTop, .after = cmBot) %>% 
  mutate(AgeBot = tl08_age_tapas$AgeBot, .after = AgeTop) %>% 
  mutate(Volume = 1, .after = AgeBot)

tl08_intercept <- median(diff(tl08_tapas$AgeTop))
tl08_pretreat <- pretreatment_data(tl08_tapas, yrInterp = tl08_intercept)

tl08_gam_pre <- tapas2mgcv(tl08_pretreat)
tl08_gam_pre <- transform(tl08_gam_pre, age = age_top - ((age_top-age_bot)/2))
tl08_gam_pre <- tl08_gam_pre[-length(tl08_gam_pre$cm_top), ] #remove the last row containing sandAR = 0
tl08_gam_pre$age <- round(tl08_gam_pre$age, digits = 3)
tl08_age_span <- max(tl08_gam_pre$age) - min(tl08_gam_pre$age)

set.seed(12)
tl08_gam <- gam(SandAR ~ s(age, k = 50), data = tl08_gam_pre, method = "REML")
tl08_gam_post <- mgcv2tapas_krzypl(tl08_gam)  
summary(tl08_gam)

tl08_thresh_gl <- global_thresh(
  series = tl08_gam_post,
  proxy = "SandAR",
  thresh.value = 0.999,
  smoothing.yr = tl08_age_span,
  keep_consecutive = FALSE,
  minCountP = 0.05,
  MinCountP_window = tl08_intercept*5
)

Plot.Anomalies(
  series = tl08_thresh_gl,
  plot.crosses = TRUE,
  plot.x = TRUE,
  plot.neg = FALSE)