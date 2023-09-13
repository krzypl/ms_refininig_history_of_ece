library(tidyverse)
library(patchwork)

wave_base_averageal <- function(F, w = 7.4, g = 9.81) {
  x <- 0.46*((g*F)/w^2)^0.28
  T <- (x*w)/g
  WB <- 0.39*T^2
  return(WB)
}

wave_base_storm <- function(F, w = 28.6, g = 9.81) {
  x <- 0.46*((g*F)/w^2)^0.28
  T <- (x*w)/g
  WB <- 0.39*T^2
  return(WB)
}

wave_base_extreme <- function(F, w = 49.6, g = 9.81) {
  x <- 0.46*((g*F)/w^2)^0.28
  T <- (x*w)/g
  WB <- 0.39*T^2
  return(WB)
}

#average refers to mean annual wind speed which approximates 7.4 m/s. Storm refers to a NOAA's threshold wind speed of storms equals to 25.9 m/s. Extreme refers to the highest wind speed historically recorded in the southern Burin, which was 47.8 m/s during the Hurricane Michael in 2000. F is a distance to the coring site.

#average---

tl27_distance2shore <- tibble(coreID = "TL27-1",
                              coring_depth = 0.5,
                              direction = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                              distance = c(46,39,118,68,51,118,83,80))
tl27_wb_average <- vector("double", length(tl27_distance2shore$distance))

for(i in seq_along(tl27_distance2shore$distance)) {
  tl27_wb_average[[i]] <- round(wave_base_averageal(tl27_distance2shore$distance[[i]]), digits = 2)
}

tl18_distance2shore <- tibble(coreID = "TL18-2",
                              coring_depth = 3.1,
                              direction = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                              distance = c(285,177,137,296,247,216,272,386))
tl18_wb_average <- vector("double", length(tl18_distance2shore$distance))

for(i in seq_along(tl18_distance2shore$distance)) {
  tl18_wb_average[[i]] <- round(wave_base_averageal(tl18_distance2shore$distance[[i]]), digits = 2)
}

tl09_distance2shore <- tibble(coreID = "TL09-1",
                              coring_depth = 1.7,
                              direction = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                              distance = c(37,116,123,39,23,30,29,37))
tl09_wb_average <- vector("double", length(tl09_distance2shore$distance))

for(i in seq_along(tl09_distance2shore$distance)) {
  tl09_wb_average[[i]] <- round(wave_base_averageal(tl09_distance2shore$distance[[i]]), digits = 2)
}

tl08_distance2shore <- tibble(coreID = "TL08-2",
                              coring_depth = 2.1,
                              direction = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                              distance = c(186,208,196,373,166,141,189,69))
tl08_wb_average <- vector("double", length(tl08_distance2shore$distance))

for(i in seq_along(tl08_distance2shore$distance)) {
  tl08_wb_average[[i]] <- round(wave_base_averageal(tl08_distance2shore$distance[[i]]), digits = 2)
}

#storm----
tl27_wb_storm <- vector("double", length(tl27_distance2shore$distance))

for(i in seq_along(tl27_distance2shore$distance)) {
  tl27_wb_storm[[i]] <- round(wave_base_storm(tl27_distance2shore$distance[[i]]), digits = 2)
}

tl18_wb_storm <- vector("double", length(tl18_distance2shore$distance))

for(i in seq_along(tl18_distance2shore$distance)) {
  tl18_wb_storm[[i]] <- round(wave_base_storm(tl18_distance2shore$distance[[i]]), digits = 2)
}

tl09_wb_storm <- vector("double", length(tl09_distance2shore$distance))

for(i in seq_along(tl09_distance2shore$distance)) {
  tl09_wb_storm[[i]] <- round(wave_base_storm(tl09_distance2shore$distance[[i]]), digits = 2)
}

tl08_wb_storm <- vector("double", length(tl08_distance2shore$distance))

for(i in seq_along(tl08_distance2shore$distance)) {
  tl08_wb_storm[[i]] <- round(wave_base_storm(tl08_distance2shore$distance[[i]]), digits = 2)
}

#extreme---
tl27_wb_extreme <- vector("double", length(tl27_distance2shore$distance))

for(i in seq_along(tl27_distance2shore$distance)) {
  tl27_wb_extreme[[i]] <- round(wave_base_extreme(tl27_distance2shore$distance[[i]]), digits = 2)
}

tl18_wb_extreme <- vector("double", length(tl18_distance2shore$distance))

for(i in seq_along(tl18_distance2shore$distance)) {
  tl18_wb_extreme[[i]] <- round(wave_base_extreme(tl18_distance2shore$distance[[i]]), digits = 2)
}

tl09_wb_extreme <- vector("double", length(tl09_distance2shore$distance))

for(i in seq_along(tl09_distance2shore$distance)) {
  tl09_wb_extreme[[i]] <- round(wave_base_extreme(tl09_distance2shore$distance[[i]]), digits = 2)
}

tl08_wb_extreme <- vector("double", length(tl08_distance2shore$distance))

for(i in seq_along(tl08_distance2shore$distance)) {
  tl08_wb_extreme[[i]] <- round(wave_base_extreme(tl08_distance2shore$distance[[i]]), digits = 2)
}

#combine records--
wb_combined <- tl27_distance2shore %>% 
  add_row(tl18_distance2shore) %>% 
  add_row(tl09_distance2shore) %>% 
  add_row(tl08_distance2shore) %>% 
  mutate(wb_average = c(tl27_wb_average,
                       tl18_wb_average,
                       tl09_wb_average,
                       tl08_wb_average),
         wb_storm = c(tl27_wb_storm,
                      tl18_wb_storm,
                      tl09_wb_storm,
                      tl08_wb_storm),
         wb_extreme = c(tl27_wb_extreme,
                      tl18_wb_extreme,
                      tl09_wb_extreme,
                      tl08_wb_extreme),
         coreID = as.factor(coreID),
         direction = as.factor(direction))

wb_combined$coreID <- factor(wb_combined$coreID, levels = c("TL27-1", "TL18-2", "TL09-1", "TL08-2"))
wb_combined$direction <- factor(wb_combined$direction, levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

wb_average_plot <- ggplot(wb_combined, aes(x = direction, y = wb_average)) +
  geom_col() +
  geom_hline(aes(yintercept = coring_depth), linewidth = 1) +
  scale_y_reverse() +
  facet_wrap(.~coreID, ncol = 1) +
  theme_bw() + 
  ylab("Depth (m)") +
  ggtitle(label = "(A) Wave bases for \n average winds")

wb_storm_plot <- ggplot(wb_combined, aes(x = direction, y = wb_storm)) +
  geom_col() +
  geom_hline(aes(yintercept = coring_depth), linewidth = 1) +
  scale_y_reverse() +
  facet_wrap(.~coreID, ncol = 1) +
  theme_bw() +
  xlab("Wind direction") +
  ggtitle(label = "(B) Wave bases for \n storm winds")

wb_extreme_plot <- ggplot(wb_combined, aes(x = direction, y = wb_extreme)) +
  geom_col() +
  geom_hline(aes(yintercept = coring_depth), linewidth = 1) +
  scale_y_reverse() +
  facet_wrap(.~coreID, ncol = 1) +
  theme_bw() +
  ggtitle(label = "(C) Wave bases for \n extreme winds") 

wb_wrapped_plot <- wrap_plots(
  wb_average_plot + 
    theme(plot.title = element_text(size = 10), axis.title.x= element_blank(), axis.text = element_text(size = 6)),
  wb_storm_plot +
    theme(plot.title = element_text(size = 10), axis.title.y = element_blank(), axis.text = element_text(size = 6)),
  wb_extreme_plot + 
    theme(plot.title = element_text(size = 10), axis.title.y = element_blank(), axis.text = element_text(size = 6), axis.title.x= element_blank()),
  ncol = 3
)

ggsave(filename="figures/wind_bases.svg", 
       plot = wb_wrapped_plot, 
       device = svg, 
       width = 130, 
       height = 120, 
       units = "mm")

ggsave(filename="figures/wind_bases.pdf", 
       plot = wb_wrapped_plot, 
       device = pdf, 
       width = 145, 
       height = 120, 
       units = "mm")

wb_diff <- wb_combined %>% 
  group_by(coreID) %>% 
  summarise(wb_average_min = min(wb_average),
            wb_average_max = max(wb_average),
            wb_extreme_min = min(wb_extreme),
            wb_extreme_max = max(wb_extreme)) %>% 
  mutate(wb_average_diff = wb_average_max - wb_average_min,
         wb_extreme_diff = wb_extreme_max - wb_extreme_min)
  