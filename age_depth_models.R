library(serac)
library(tidyverse)


#########################################TL27-1 core####################
TL27 <- serac(name="TL27",
              model = c("CFCS", "CIC", "CRS"),
              coring_yr = 2019,
              plot_Pb=TRUE,
              plot_Cs=TRUE,
              plotpdf = TRUE,
              Hemisphere=c("NH"),
              FF=c(160,200),
              mass_depth = TRUE)

######extrapolation of ages based on mass depth
tl27_cfcs <- TL27$`CFCS age-depth model interpolated` #extract interpolated ages (CFCS model was used as the most appropriate)
tl27_sed <- read_csv("data/tl27_sed.csv") #read dataset for TL27

#prepare object for extrapolation
tl27_age_ad_prep <- approx(tl27_cfcs$depth_avg_mm/10, tl27_cfcs$BestAD,
                           xout = tl27_sed$depth) 
tl27_age_ad_prep_unique <- unique(tl27_age_ad_prep$y) #get rid of replicated dates at the end of the age column
tl27_age_ad <- tibble(ageAD = c(tl27_age_ad_prep_unique,
                                rep(NA, times = length(tl27_age_ad_prep$x)-length(tl27_age_ad_prep_unique))),
                      depth = tl27_age_ad_prep$x) 

tl27_mass_depth <- c(unique(tl27_cfcs$mass_depth_g.cm.2[-1]), #get rid of replicated mass depth values at the end of the column
                     rep(NA, times = (length(tl27_sed$DBD) - length(unique(tl27_cfcs$mass_depth_g.cm.2[-1]))))) 

tl27_age_ad <- tl27_age_ad %>% 
  mutate(DBD = tl27_sed$DBD, mass_depth = tl27_mass_depth)

tl27_sel_dbd <- tl27_age_ad$DBD[is.na(tl27_age_ad$mass_depth)] * 0.5 #DBD values are multiplied by 0.5 for further calculation because sampling resolution was 0.5 cm
tl27_sel_mass_depth <- tl27_age_ad$mass_depth[is.na(tl27_age_ad$mass_depth)]
tl27_sel_mass_depth[1] <- max(tl27_age_ad$mass_depth, na.rm = TRUE) + tl27_sel_dbd[1]

#creating a loop for calculating all the missing mass depth values
tl27_out_mass_depth <- vector("double", length(tl27_sel_mass_depth)) 
tl27_out_mass_depth[1] <- tl27_sel_mass_depth[1]
for (i in 2:length(tl27_sel_mass_depth)) {
  tl27_out_mass_depth[[i]] <- tl27_out_mass_depth[[i-1]] + tl27_sel_dbd[[i]]
}

tl27_age_ad$mass_depth[is.na(tl27_age_ad$mass_depth)] <- tl27_out_mass_depth

#performing the extrapolation with the calculated mass depth values
tl27_age_ad$ageAD[is.na(tl27_age_ad$ageAD)] <- 2019 - tl27_age_ad$mass_depth[is.na(tl27_age_ad$ageAD)]/abs(TL27$`CFCS mass accumulation rate`$MAR_g.cm.2.yr.1)

write_csv(tl27_age_ad, "data/tl27_age_ad.csv")

######age estimates for peak detection analysis
tl27_age_tapas_top <- approx(tl27_age_ad$depth, tl27_age_ad$ageAD,
                           xout = tl27_sed$cmTop)$y
tl27_age_tapas_top[1] <- 2019 #age for 0 cm
tl27_age_tapas_top <- 1950 - tl27_age_tapas_top

tl27_age_tapas_bot <- approx(tl27_age_ad$depth, tl27_age_ad$ageAD,
                             xout = tl27_sed$cmBot)$y

tl27_age_tapas_bot[max(length(tl27_age_tapas_bot))] <- tl27_age_tapas_bot[max(length(tl27_age_tapas_bot)) - 1] - (tl27_age_tapas_bot[max(length(tl27_age_tapas_bot)) - 2] - tl27_age_tapas_bot[max(length(tl27_age_tapas_bot)) - 1]) #calculate approximate age of the bottom of the last sample

tl27_age_tapas_bot <- 1950 - tl27_age_tapas_bot

tl27_age_tapas <- tibble(AgeTop = tl27_age_tapas_top, AgeBot = tl27_age_tapas_bot)
write_csv(tl27_age_tapas, "data/tl27_age_tapas.csv")


#########################################TL18-2 core#####################
TL18 <- serac(name="TL18",
              model = c("CFCS","CRS","CIC"),
              coring_yr = 2019,
              Hemisphere= c("NH"),
              plot_Cs = TRUE,
              plot_Pb = TRUE,
              plotpdf = TRUE,
              FF=c(100,140),
              NWT=c(80,100),
              mass_depth = TRUE)

######extrapolation of ages based on mass depth
tl18_cfcs <- TL18$`CFCS age-depth model interpolated` #extract interpolated ages (CFCS model was used as the most appropriate)
tl18_sed <- read_csv("data/tl18_sed.csv") #read dataset for TL18

#prepare object for extrapolation
tl18_age_ad_prep <- approx(tl18_cfcs$depth_avg_mm/10, tl18_cfcs$BestAD,
                           xout = tl18_sed$depth) 
tl18_age_ad_prep_unique <- unique(tl18_age_ad_prep$y) #get rid of replicated dates at the end of the age column
tl18_age_ad <- tibble(ageAD = c(tl18_age_ad_prep_unique,
                                rep(NA, times = length(tl18_age_ad_prep$x)- length(tl18_age_ad_prep_unique))),
                      depth = tl18_age_ad_prep$x) 

tl18_mass_depth <- c(unique(tl18_cfcs$mass_depth_g.cm.2[-1]), #get rid of replicated mass depth values at the end of the column
                     rep(NA, times = (length(tl18_sed$DBD) - length(unique(tl18_cfcs$mass_depth_g.cm.2[-1]))))) 

tl18_age_ad <- tl18_age_ad %>% 
  mutate(DBD = tl18_sed$DBD, mass_depth = tl18_mass_depth)

tl18_sel_dbd <- tl18_age_ad$DBD[is.na(tl18_age_ad$mass_depth)] * 0.5 #DBD values are multiplied by 0.5 for further calculation because sampling resolution was 0.5 cm
tl18_sel_mass_depth <- tl18_age_ad$mass_depth[is.na(tl18_age_ad$mass_depth)]
tl18_sel_mass_depth[1] <- max(tl18_age_ad$mass_depth, na.rm = TRUE) + tl18_sel_dbd[1]

#creating a loop for calculating all the missing mass depth values
tl18_out_mass_depth <- vector("double", length(tl18_sel_mass_depth)) 
tl18_out_mass_depth[1] <- tl18_sel_mass_depth[1]
for (i in 2:length(tl18_sel_mass_depth)) {
  tl18_out_mass_depth[[i]] <- tl18_out_mass_depth[[i-1]] + tl18_sel_dbd[[i]]
}

tl18_age_ad$mass_depth[is.na(tl18_age_ad$mass_depth)] <- tl18_out_mass_depth

#performing the extrapolation with the calculated mass depth values
tl18_age_ad$ageAD[is.na(tl18_age_ad$ageAD)] <- 2019 - tl18_age_ad$mass_depth[is.na(tl18_age_ad$ageAD)]/abs(TL18$`CFCS mass accumulation rate`$MAR_g.cm.2.yr.1)

write_csv(tl18_age_ad, "data/tl18_age_ad.csv")

######age estimates for peak detection analysis
tl18_age_tapas_top <- approx(tl18_age_ad$depth, tl18_age_ad$ageAD,
                             xout = tl18_sed$cmTop)$y
tl18_age_tapas_top[1] <- 2019 #age for 0 cm
tl18_age_tapas_top <- 1950 - tl18_age_tapas_top

tl18_age_tapas_bot <- approx(tl18_age_ad$depth, tl18_age_ad$ageAD,
                             xout = tl18_sed$cmBot)$y

tl18_age_tapas_bot[max(length(tl18_age_tapas_bot))] <- tl18_age_tapas_bot[max(length(tl18_age_tapas_bot)) - 1] - (tl18_age_tapas_bot[max(length(tl18_age_tapas_bot)) - 2] - tl18_age_tapas_bot[max(length(tl18_age_tapas_bot)) - 1]) #calculate approximate age of the bottom of the last sample

tl18_age_tapas_bot <- 1950 - tl18_age_tapas_bot

tl18_age_tapas <- tibble(AgeTop = tl18_age_tapas_top, AgeBot = tl18_age_tapas_bot)
write_csv(tl18_age_tapas, "data/tl18_age_tapas.csv")

#########################################TL09-1 core##################
TL09 <- serac(name="TL09",
              model = c("CFCS","CRS","CIC","CRS_pw"),
              coring_yr = 2019,
              plot_Pb=TRUE,
              plot_Cs=TRUE,
              plotpdf = TRUE,
              Hemisphere= "NH",
              FF=c(220,260),
              NWT=c(180,220),
              age_forced_CRS=c(1929),
              depth_forced_CRS =c(310),
              #inst_deposit = c(311,330), 
              ignore = 325,
              mass_depth = TRUE) #(the figure in the manuscript was produced with inst_depost = c(311,330)#

######extrapolation of ages based on mass depth
tl09_crspw <- TL09$`CRS piecewise model age-depth model interpolated` #extract interpolated ages (CRS piecewise model was used as the most appropriate)
tl09_sed <- read_csv("data/tl09_sed.csv") #read dataset for TL09

#prepare object for extrapolation
tl09_age_ad_prep <- approx(tl09_crspw$depth_avg_mm/10, tl09_crspw$BestAD,
                           xout = tl09_sed$depth) 
tl09_age_ad <- tibble(ageAD = tl09_age_ad_prep$y, depth = tl09_age_ad_prep$x)
tl09_age_ad$ageAD[is.na(tl09_age_ad$ageAD)] <- 1929

write_csv(tl09_age_ad, "data/tl09_age_ad.csv")

######age estimates for peak detection analysis
tl09_age_tapas_top <- approx(tl09_age_ad$depth, tl09_age_ad$ageAD,
                             xout = tl09_sed$cmTop)$y
tl09_age_tapas_top[1] <- 2019 #age for 0 cm
tl09_age_tapas_top <- 1950 - tl09_age_tapas_top

tl09_age_tapas_bot <- approx(tl09_age_ad$depth, tl09_age_ad$ageAD,
                             xout = tl09_sed$cmBot)$y

tl09_age_tapas_bot[max(length(tl09_age_tapas_bot))] <- tl09_age_tapas_bot[max(length(tl09_age_tapas_bot)) - 1] - (tl09_age_tapas_bot[max(length(tl09_age_tapas_bot)) - 2] - tl09_age_tapas_bot[max(length(tl09_age_tapas_bot)) - 1]) #calculate approximate age of the bottom of the last sample

tl09_age_tapas_bot <- 1950 - tl09_age_tapas_bot

tl09_age_tapas <- tibble(AgeTop = tl09_age_tapas_top, AgeBot = tl09_age_tapas_bot)
write_csv(tl09_age_tapas, "data/tl09_age_tapas.csv")

#########################################TL08-2 core###################
TL08 <- serac(name="TL08",
              model = c("CFCS","CRS","CIC"),
              coring_yr = 2019,
              plot_Pb=TRUE,
              plot_Cs=TRUE,
              plotpdf = TRUE,
              Hemisphere=c("NH"),
              FF=c(200,240),
              NWT=c(160,190),
              mass_depth = TRUE)

######extrapolation of ages based on mass depth
tl08_cfcs <- TL08$`CFCS age-depth model interpolated` #extract interpolated ages (CFCS model was used as the most appropriate)
tl08_sed <- read_csv("data/tl08_sed.csv") #read dataset for TL08

#prepare object for extrapolation
tl08_age_ad_prep <- approx(tl08_cfcs$depth_avg_mm/10, tl08_cfcs$BestAD,
                           xout = tl08_sed$depth) 
tl08_age_ad_prep_unique <- unique(tl08_age_ad_prep$y) #get rid of replicated dates at the end of the age column
tl08_age_ad <- tibble(ageAD = c(tl08_age_ad_prep_unique,
                                rep(NA, times = length(tl08_age_ad_prep$x)- length(tl08_age_ad_prep_unique))),
                      depth = tl08_age_ad_prep$x) 

tl08_mass_depth <- c(unique(tl08_cfcs$mass_depth_g.cm.2[-1]), #get rid of replicated mass depth values at the end of the column
                     rep(NA, times = (length(tl08_sed$DBD) - length(unique(tl08_cfcs$mass_depth_g.cm.2[-1]))))) 

tl08_age_ad <- tl08_age_ad %>% 
  mutate(DBD = tl08_sed$DBD, mass_depth = tl08_mass_depth)

tl08_sel_dbd <- tl08_age_ad$DBD[is.na(tl08_age_ad$mass_depth)] * 0.5 #DBD values are multiplied by 0.5 for further calculation because sampling resolution was 0.5 cm
tl08_sel_mass_depth <- tl08_age_ad$mass_depth[is.na(tl08_age_ad$mass_depth)]
tl08_sel_mass_depth[1] <- max(tl08_age_ad$mass_depth, na.rm = TRUE) + tl08_sel_dbd[1]

#creating a loop for calculating all the missing mass depth values
tl08_out_mass_depth <- vector("double", length(tl08_sel_mass_depth)) 
tl08_out_mass_depth[1] <- tl08_sel_mass_depth[1]
for (i in 2:length(tl08_sel_mass_depth)) {
  tl08_out_mass_depth[[i]] <- tl08_out_mass_depth[[i-1]] + tl08_sel_dbd[[i]]
}

tl08_age_ad$mass_depth[is.na(tl08_age_ad$mass_depth)] <- tl08_out_mass_depth

#performing the extrapolation with the calculated mass depth values
tl08_age_ad$ageAD[is.na(tl08_age_ad$ageAD)] <- 2019 - tl08_age_ad$mass_depth[is.na(tl08_age_ad$ageAD)]/abs(TL08$`CFCS mass accumulation rate`$MAR_g.cm.2.yr.1)

write_csv(tl08_age_ad, "data/tl08_age_ad.csv")

######age estimates for peak detection analysis
tl08_age_tapas_top <- approx(tl08_age_ad$depth, tl08_age_ad$ageAD,
                             xout = tl08_sed$cmTop)$y
tl08_age_tapas_top[1] <- 2019 #age for 0 cm
tl08_age_tapas_top <- 1950 - tl08_age_tapas_top

tl08_age_tapas_bot <- approx(tl08_age_ad$depth, tl08_age_ad$ageAD,
                             xout = tl08_sed$cmBot)$y

tl08_age_tapas_bot[max(length(tl08_age_tapas_bot))] <- tl08_age_tapas_bot[max(length(tl08_age_tapas_bot)) - 1] - (tl08_age_tapas_bot[max(length(tl08_age_tapas_bot)) - 2] - tl08_age_tapas_bot[max(length(tl08_age_tapas_bot)) - 1]) #calculate approximate age of the bottom of the last sample

tl08_age_tapas_bot <- 1950 - tl08_age_tapas_bot

tl08_age_tapas <- tibble(AgeTop = tl08_age_tapas_top, AgeBot = tl08_age_tapas_bot)
write_csv(tl08_age_tapas, "data/tl08_age_tapas.csv")


#for the manuscript, each model was rerun to obtain 137Cs profile in depth scale with mass_depth argument set to FALSE