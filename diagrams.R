library(tidyverse)
library(tidypaleo)
library(patchwork)
library(rioja)
library(vegan)
library(ggvegan)
source("diatom_eco_groups.R")
theme_set(theme_paleo(base_size = 9))
          


###########tl27 sed + diatom###############
tl27_sed <- read_csv("data/tl27_sed.csv")
tl27_age_ad <- read_csv("data/tl27_age_ad.csv")

tl27_sed <- tl27_sed %>% 
  mutate(ageAD = tl27_age_ad$ageAD, .after = location)

tl27_sed <- tl27_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl27_tidy <- tl27_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, ageAD, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl27_ages <- tl27_sed %>% 
  select(depth, ageAD)

tl27_diatom <- read_csv("data/tl27_diatom.csv")

tl27_diatom_depth <- tl27_diatom$depth
tl27_diatom_spec <- tl27_diatom[,-1]
tl27_diatom_perc <- tl27_diatom_spec/rowSums(tl27_diatom_spec, na.rm = TRUE)*100

tl27_diatom_eco <- diatom_eco_groups(tl27_diatom_perc, tl27_diatom_depth) %>% 
  add_row(tibble(depth = max(tl27_sed$depth), M = 0, B = 0, E = 0, F = 0, U = 0))

tl27_diatom_eco_tidy <- tl27_diatom_eco %>% 
  pivot_longer(!depth, names_to = "ecological_group", values_to = "rel_abund") %>% 
  mutate(ecological_group = fct_relevel(ecological_group, "M", "B", "E", "F", "U"))

tl27_plot <- ggplot(tl27_tidy, aes(x = value, y = depth)) +
  geom_lineh(linewidth = 0.05) +
  geom_point(size = 0.05) +
  scale_y_reverse(breaks = seq(0, 35, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(A) TL27, TL27-1")

tl27_adm <- age_depth_model(
  tl27_ages,
  depth = depth,
  age = ageAD
)

tl27_diatom_eco_tidy_ex <- tl27_diatom_eco_tidy %>% 
  filter(ecological_group %in% c("M", "B", "E", "U")) %>% 
  mutate(rel_abund = rel_abund*5)

tl27_diatom_plot <- ggplot(tl27_diatom_eco_tidy, aes(x = rel_abund, y = depth)) +
  geom_col_segsh() +
  geom_col_segsh(data = filter(tl27_diatom_eco_tidy, depth %in% c(20.75, 35.25)),
                               aes(x = rel_abund, y = depth),
                               color = "red") +
  scale_y_reverse(breaks = seq(0, 35, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_abundanceh(vars(ecological_group), rotate_facet_labels = 0, dont_italicize = c("M", "B", "E", "F", "U")) +
  geom_lineh(data = tl27_diatom_eco_tidy_ex, aes(x = rel_abund, y = depth), col = "grey70", lty = 2, lwd = 0.2) +
  labs(subtitle = "Diatom ecological groups (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))

tl27_wrapped <- wrap_plots(
  tl27_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()) +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)),
  tl27_diatom_plot +
    scale_y_depth_age(
      tl27_adm,
      age_name = "CFCS-derived age (AD)",
      breaks = c(breaks = seq(0, 35, by = 5)),
      expand = expansion(mult = c(0.02, 0.02)),
      age_breaks = c(2019, seq(1920, 2000, by = 20)),
      age_labels = as.character(c(2019, seq(1920, 2000, by = 20)))
    ) +
    theme(axis.text.y.left = element_blank(), axis.title.x.bottom =  element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 3)
)  

###########tl18_2 sed + diatom##########
tl18_sed <- read_csv("data/tl18_sed.csv")
tl18_age_ad <- read_csv("data/tl18_age_ad.csv")
tl18_sed <- tl18_sed %>% 
  mutate(ageAD = tl18_age_ad$ageAD, .after = location)

tl18_sed <- tl18_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl18_tidy <- tl18_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, ageAD, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl18_ages <- tl18_sed %>% 
  select(depth, ageAD)

tl18_diatom <- read_csv("data/tl18_diatom.csv")

tl18_diatom_depth <- tl18_diatom$depth
tl18_diatom_spec <- tl18_diatom[,-1]
tl18_diatom_perc <- tl18_diatom_spec/rowSums(tl18_diatom_spec, na.rm = TRUE)*100

tl18_diatom_eco <- diatom_eco_groups(tl18_diatom_perc, tl18_diatom_depth) %>% 
  add_row(tibble(depth = max(tl18_sed$depth), M = 0, B = 0, E = 0, F = 0, U = 0))

tl18_diatom_eco_tidy <- tl18_diatom_eco %>% 
  pivot_longer(!depth, names_to = "ecological_group", values_to = "rel_abund") %>% 
  mutate(ecological_group = fct_relevel(ecological_group, "M", "B", "E", "F", "U"))

tl18_plot <- ggplot(tl18_tidy, aes(x = value, y = depth)) +
  geom_lineh(size = 0.05) +
  geom_point(size = 0.05) +
  geom_point(data = filter(tl18_tidy, depth %in% c(5.75, 23.75, 31.75, 35.75, 40.25) & param == 'Sand >0.25 mm'),
             shape = 8, color = "red", size = 2) +
  scale_y_reverse(breaks = seq(0, 40, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(B) Broad Pond, TL18-2")

tl18_adm <- age_depth_model(
  tl18_ages,
  depth = depth,
  age = ageAD
)

tl18_diatom_eco_tidy_ex <- tl18_diatom_eco_tidy %>% 
  filter(ecological_group %in% c("M", "B", "E", "U")) %>% 
  mutate(rel_abund = rel_abund*5)

tl18_diatom_plot <- ggplot(tl18_diatom_eco_tidy, aes(x = rel_abund, y = depth)) +
  geom_col_segsh() +
  geom_col_segsh(data = filter(tl18_diatom_eco_tidy, depth %in% c(5.75, 23.75, 31.75, 35.75, 40.25)),
                 aes(x = rel_abund, y = depth),
                 color = "red") +
  scale_y_reverse(breaks = seq(0, 40, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_abundanceh(vars(ecological_group), rotate_facet_labels = 0, dont_italicize = c("M", "B", "E", "F", "U")) +
  geom_lineh(data = tl18_diatom_eco_tidy_ex, aes(x = rel_abund, y = depth), col = "grey70", lty = 2, lwd = 0.2) +
  labs(subtitle = "Diatom ecological groups (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))

tl18_2_wrapped <- wrap_plots(
  tl18_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()) +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)
    ),
  tl18_diatom_plot +
    scale_y_depth_age(
      tl18_adm,
      age_name = "CFCS-derived age (AD)",
      breaks = c(breaks = seq(0, 40, by = 5)),
      expand = expansion(mult = c(0.02, 0.02)),
      age_breaks = c(2019, seq(1640, 2000, by = 40)),
      age_labels = as.character(c(2019, seq(1640, 2000, by = 40)))
    ) +
    theme(axis.text.y.left = element_blank(), axis.title.x.bottom =  element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 3)
)

###########tl18_1 sed####
tl18_1_sed <- read_csv("data/tl18_1_sed.csv")

tl18_1_sed <- tl18_1_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl18_1_tidy <- tl18_1_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl18_1_plot <- ggplot(tl18_1_tidy, aes(x = value, y = depth)) +
  geom_lineh(size = 0.05) +
  geom_point(size = 0.05) +
  scale_y_reverse(breaks = seq(0, 40, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(B) Broad Pond, TL18-1")


###########tl09_1 sed + diatom##########
tl09_sed <- read_csv("data/tl09_sed.csv")
tl09_age_ad <- read_csv("data/tl09_age_ad.csv")
tl09_sed <- tl09_sed %>% 
  mutate(ageAD = tl09_age_ad$ageAD, .after = location)

tl09_sed <- tl09_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl09_tidy <- tl09_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, ageAD, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl09_ages <- tl09_sed %>% 
  select(depth, ageAD)

#creating of age axis with tidypaleo requires values to change monotonically, so the two bottommost samples have been modified. Note, the change does not influence the position of tick marks on the age scale
tl09_ages$ageAD[32] <- 1928
tl09_ages$ageAD[33] <- 1927

tl09_diatom <- read_csv("data/tl09_diatom.csv")

tl09_diatom_depth <- tl09_diatom$depth
tl09_diatom_spec <- tl09_diatom[,-1]
tl09_diatom_perc <- tl09_diatom_spec/rowSums(tl09_diatom_spec, na.rm = TRUE)*100

tl09_diatom_eco <- diatom_eco_groups(tl09_diatom_perc, tl09_diatom_depth)

tl09_diatom_eco_tidy <- tl09_diatom_eco %>% 
  pivot_longer(!depth, names_to = "ecological_group", values_to = "rel_abund") %>% 
  mutate(ecological_group = fct_relevel(ecological_group, "M", "B", "E", "F", "U"))
  
tl09_plot <- ggplot(tl09_tidy, aes(x = value, y = depth)) +
  geom_lineh(linewidth = 0.05) +
  geom_point(size = 0.05) +
  geom_point(data = filter(tl09_tidy, depth == 27.5 & param == 'Sand >0.25 mm'), shape = 8, color = "red", size = 2) +
  scale_y_reverse(breaks = seq(0, 30, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(C) Porsh Pond, TL09-1")

tl09_adm <- age_depth_model(
  tl09_ages,
  depth = depth,
  age = ageAD
)

tl09_diatom_eco_tidy_ex <- tl09_diatom_eco_tidy %>% 
  filter(ecological_group %in% c("M", "B", "E", "U")) %>% 
  mutate(rel_abund = rel_abund*5)

tl09_diatom_plot <- ggplot(tl09_diatom_eco_tidy, aes(x = rel_abund, y = depth)) +
  geom_col_segsh() +
  scale_y_reverse(breaks = seq(0, 30, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_abundanceh(vars(ecological_group), rotate_facet_labels = 0,
                   dont_italicize = c("M", "B", "E", "F", "U"))+
  geom_col_segsh(data = filter(tl09_diatom_eco_tidy, depth %in% c(27.5, 31.5, 32.5)),
                 aes(x = rel_abund, y = depth), color = "red") +
  geom_lineh(data = tl09_diatom_eco_tidy_ex, aes(x = rel_abund, y = depth), col = "grey70", lty = 2, lwd = 0.2) +
  labs(subtitle = "Diatom ecological groups (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))
  

tl09_wrapped <- wrap_plots(
  tl09_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()) +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)
    ),
  tl09_diatom_plot +
    scale_y_depth_age(
      tl09_adm,
      age_name = "CRS piecewise-derived age (AD)",
      breaks = c(breaks = seq(0, 30, by = 5)),
      expand = expansion(mult = c(0.02, 0.02)),
      age_breaks = c(2019, seq(1930, 2000, by = 10)),
      age_labels = as.character(c(2019, seq(1930, 2000, by = 10)))
    ) +
    theme(axis.text.y.left = element_blank(), axis.title.x.bottom =  element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 3)
)

#######tl08_2 sed + diatom########
tl08_sed <- read_csv("data/tl08_sed.csv")
tl08_age_ad <- read_csv("data/tl08_age_ad.csv")

tl08_sed <- tl08_sed %>% 
  mutate(ageAD = tl08_age_ad$ageAD, .after = location)

tl08_sed <- tl08_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl08_tidy <- tl08_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, ageAD, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl08_ages <- tl08_sed %>% 
  select(depth, ageAD)

tl08_diatom <- read_csv("data/tl08_diatom.csv")
tl08_diatom_depth <- tl08_diatom$depth
tl08_diatom_spec <- tl08_diatom[,-1]
tl08_diatom_perc <- tl08_diatom_spec/rowSums(tl08_diatom_spec, na.rm = TRUE)*100

tl08_diatom_eco <- diatom_eco_groups(tl08_diatom_perc, tl08_diatom_depth) %>% 
  add_row(tibble(depth = max(tl08_sed$depth), M = 0, B = 0, E = 0, F = 0, U = 0))

tl08_diatom_eco_tidy <- tl08_diatom_eco %>% 
  pivot_longer(!depth, names_to = "ecological_group", values_to = "rel_abund") %>% 
  mutate(ecological_group = fct_relevel(ecological_group, "M", "B", "E", "F", "U"))

tl08_plot <- ggplot(tl08_tidy, aes(x = value, y = depth)) +
  geom_lineh(linewidth = 0.05) +
  geom_point(size = 0.05) +
  geom_point(data = filter(tl08_tidy, depth == 27.75 & param == 'Sand >0.25 mm'), shape = 8, color = "red", size = 2) +
  scale_y_reverse(breaks = seq(0, 30, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(C) Frenchmans Pond, TL08-2")

tl08_adm <- age_depth_model(
  tl08_ages,
  depth = depth,
  age = ageAD
)

tl08_diatom_eco_tidy_ex <- tl08_diatom_eco_tidy %>% 
  filter(ecological_group %in% c("M", "B", "E", "U")) %>% 
  mutate(rel_abund = rel_abund*5)

tl08_diatom_plot <- ggplot(tl08_diatom_eco_tidy, aes(x = rel_abund, y = depth)) +
  geom_col_segsh() +
  scale_y_reverse(breaks = seq(0, 30, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_abundanceh(vars(ecological_group), rotate_facet_labels = 0, dont_italicize = c("M", "B", "E", "F", "U")) +
  geom_col_segsh(data = filter(tl08_diatom_eco_tidy, depth == 27.75), aes(x = rel_abund, y = depth), color = "red") +
  geom_lineh(data = tl08_diatom_eco_tidy_ex, aes(x = rel_abund, y = depth), col = "grey70", lty = 2, lwd = 0.2) +
  labs(subtitle = "Diatom ecological groups (%)") +
  theme(plot.subtitle = element_text(hjust = 0.5))

tl08_2_wrapped <- wrap_plots(
  tl08_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()) +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)
    ),
  tl08_diatom_plot +
    scale_y_depth_age(
      tl08_adm,
      age_name = "CFCS-derived age (AD)",
      breaks = c(breaks = seq(0, 30, by = 5)),
      expand = expansion(mult = c(0.02, 0.02)),
      age_breaks = c(2019, seq(1900, 2010, by = 20)),
      age_labels = as.character(c(2019, seq(1900, 2010, by = 20)))
    ) +
    theme(axis.text.y.left = element_blank(), axis.title.x.bottom =  element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 3)
)

###########tl08_1 sed####
tl08_1_sed <- read_csv("data/tl08_1_sed.csv")

tl08_1_sed <- tl08_1_sed %>% 
  rename("Sand >0.25 mm" = "Sand")

tl08_1_tidy <- tl08_1_sed %>%
  pivot_longer(cols = DBD:"Sand >0.25 mm", names_to = "param", values_to = "value") %>% 
  relocate(location, param, depth, value) %>% 
  arrange(param) %>% 
  add_column(units = NA) %>% 
  mutate(units = replace(units, is.na(units) & param == "DBD", "g cm-3")) %>% 
  mutate(units = replace(units, is.na(units) & param == "LOI", "%")) %>% 
  mutate(units = replace(units, is.na(units) & param == "Sand >0.25 mm", "count"))

tl08_1_plot <- ggplot(tl08_1_tidy, aes(x = value, y = depth)) +
  geom_lineh(size = 0.05) +
  geom_point(size = 0.05) +
  scale_y_reverse(breaks = seq(0, 40, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  ggtitle(label = "(B) Frenchmans Pond, TL08-1")


###########wrapped diagrams##########
diagrams_wrapped <- wrap_plots(
  tl27_wrapped,
  tl18_2_wrapped,
  tl09_wrapped,
  tl08_2_wrapped,
  nrow = 4
)

ggsave(filename="figures/fig4_diagrams_wrapped.svg", 
       plot = diagrams_wrapped, 
       device = svg, 
       width = 130, 
       height = 237, 
       units = "mm")

supp_sed_diagrams_wrapped <- wrap_plots(
  tl18_1_plot +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)
    ),
  tl08_1_plot +
    facet_geochem_gridh(
      vars(param),
      units = c("DBD" = "g cm-3", "LOI" = "%", "Sand >0.25 mm" = NA)
    )
)

ggsave(filename="figures/fig5_supp_sed_diagrams.svg", 
       plot = supp_sed_diagrams_wrapped, 
       device = svg, 
       width = 130, 
       height = 70, 
       units = "mm")

#######tl27 diatom supplement#########################

tl27_diatom_supp <- tl27_diatom
colnames(tl27_diatom_supp)<-gsub("[(M)]","",colnames(tl27_diatom_supp))
colnames(tl27_diatom_supp)<-gsub("B ","",colnames(tl27_diatom_supp))
colnames(tl27_diatom_supp)<-gsub("E ","",colnames(tl27_diatom_supp))
colnames(tl27_diatom_supp)<-gsub("F ","",colnames(tl27_diatom_supp))
colnames(tl27_diatom_supp)<-gsub("I ","",colnames(tl27_diatom_supp))

tl27_depth <- tl27_diatom_supp$depth
tl27_spec <- tl27_diatom_supp[,-1]
tl27_perc <- tl27_spec/rowSums(tl27_spec, na.rm = TRUE)*100
tl27_perc[is.na(tl27_perc)] <- 0

tl27_colsums <- colSums(tl27_perc)
tl27_sorted <- tl27_colsums[order(tl27_colsums, decreasing = TRUE)]
tl27_35_prep <- tl27_sorted[1:35]
tl27_nm <- names(tl27_35_prep)
tl27_pattern <- str_c(tl27_nm, collapse = "|")
tl27_35 <- select(tl27_perc, matches(tl27_pattern))

svg("figures/tl27_diatom_diagram.svg",
    width = 12.5,
    height = 5,
    bg = "white")
tl27_stratplot <- strat.plot(tl27_35,
                             title = "(A) TL27, TL27-1",
                             cex.title = 1,
                             cex.ylabel = 0.8,
                             ylabel = "Depth (cm)",
                             yvar = tl27_depth,
                             cex.xlabel = 0.8,
                             xRight = 0.95,
                             yTop = 0.6,
                             yBottom = 0.1,
                             srt.xlabel = 45,
                             y.tks = seq(0, 36, by = 2),
                             scale.percent = TRUE,
                             scale.minmax = TRUE,
                             plot.bar = TRUE,
                             plot.poly = TRUE,
                             plot.line = FALSE,
                             lwd.bar = 5,
                             col.bar = c(rep("darkred", 1),
                                         rep("darkblue", 6),
                                         rep("darkgray", 5),
                                         rep("darkgreen", 22),
                                         rep("black", 1)),
                             ylim = c(-1,36),
                             y.rev = TRUE)
addZone(tl27_stratplot, upper = c(20.75, 35.25), lty = 1, lwd = 4, col = alpha("red", 0.2))
dev.off()

#######tl18 diatom supplement#########################
tl18_diatom_supp <- tl18_diatom
colnames(tl18_diatom_supp)<-gsub("[(M)]","",colnames(tl18_diatom_supp))
colnames(tl18_diatom_supp)<-gsub("B ","",colnames(tl18_diatom_supp))
colnames(tl18_diatom_supp)<-gsub("E ","",colnames(tl18_diatom_supp))
colnames(tl18_diatom_supp)<-gsub("F ","",colnames(tl18_diatom_supp))
colnames(tl18_diatom_supp)<-gsub("I ","",colnames(tl18_diatom_supp))

tl18_depth <- tl18_diatom_supp$depth
tl18_spec <- tl18_diatom_supp[,-1]
tl18_perc <- tl18_spec/rowSums(tl18_spec, na.rm = TRUE)*100
tl18_perc[is.na(tl18_perc)] <- 0

tl18_colsums <- colSums(tl18_perc)
tl18_sorted <- tl18_colsums[order(tl18_colsums, decreasing = TRUE)]
tl18_35_prep <- tl18_sorted[1:35]
tl18_nm <- names(tl18_35_prep)
tl18_pattern <- str_c(tl18_nm, collapse = "|")
tl18_35 <- select(tl18_perc, matches(tl18_pattern))

svg("figures/tl18_diatom_diagram.svg",
    width = 12.5,
    height = 5,
    bg = "white")
tl18_stratplot <- strat.plot(tl18_35,
                             title = "(B) Broad Pond, TL18-2",
                             cex.title = 1,
                             cex.ylabel = 0.8,
                             ylabel = "Depth (cm)",
                             yvar = tl18_depth,
                             cex.xlabel = 0.8,
                             xRight = 0.95,
                             yTop = 0.6,
                             yBottom = 0.1,
                             srt.xlabel = 45,
                             y.tks = seq(0, 42, by = 2),
                             scale.percent = TRUE,
                             scale.minmax = TRUE,
                             plot.bar = TRUE,
                             plot.poly = TRUE,
                             plot.line = FALSE,
                             lwd.bar = 3,
                             col.bar = c(rep("darkred", 1),
                                         rep("darkblue", 5),
                                         rep("darkgray", 4),
                                         rep("darkgreen", 24),
                                         rep("black", 1)),
                             ylim = c(-1,43),
                             y.rev = TRUE)
addZone(tl18_stratplot, upper = c(5.75, 23.75, 31.75, 35.75, 40.25), lty = 1, lwd = 2, col = alpha("red", 0.2))
dev.off()

#######tl09 diatom supplement#########################
tl09_diatom_supp <- tl09_diatom
colnames(tl09_diatom_supp)<-gsub("[(M)]","",colnames(tl09_diatom_supp))
colnames(tl09_diatom_supp)<-gsub("B ","",colnames(tl09_diatom_supp))
colnames(tl09_diatom_supp)<-gsub("E ","",colnames(tl09_diatom_supp))
colnames(tl09_diatom_supp)<-gsub("F ","",colnames(tl09_diatom_supp))
colnames(tl09_diatom_supp)<-gsub("I ","",colnames(tl09_diatom_supp))

tl09_depth <- tl09_diatom_supp$depth
tl09_spec <- tl09_diatom_supp[,-1]
tl09_perc <- tl09_spec/rowSums(tl09_spec, na.rm = TRUE)*100
tl09_perc[is.na(tl09_perc)] <- 0

tl09_colsums <- colSums(tl09_perc)
tl09_sorted <- tl09_colsums[order(tl09_colsums, decreasing = TRUE)]
tl09_35_prep <- tl09_sorted[1:35]
tl09_nm <- names(tl09_35_prep)
tl09_pattern <- str_c(tl09_nm, collapse = "|")
tl09_35 <- select(tl09_perc, matches(tl09_pattern))

svg("figures/tl09_diatom_diagram.svg",
    width = 12.5,
    height = 5,
    bg = "white")
tl09_stratplot <- strat.plot(tl09_35,
                             cex.ylabel = 0.8,
                             title = "(C) Porsh Pond, TL09-1",
                             cex.title = 1,
                             ylabel = "Depth (cm)",
                             yvar = tl09_depth,
                             cex.xlabel = 0.8,
                             xRight = 0.95,
                             yTop = 0.6,
                             yBottom = 0.1,
                             srt.xlabel = 45,
                             y.tks = seq(0, 36, by = 2),
                             scale.percent = TRUE,
                             scale.minmax = TRUE,
                             plot.bar = TRUE,
                             plot.poly = TRUE,
                             plot.line = FALSE,
                             lwd.bar = 5,
                             col.bar = c(rep("darkred", 1),
                                         rep("darkblue", 4),
                                         rep("darkgray", 1),
                                         rep("darkgreen", 28),
                                         rep("black", 1)),
                             ylim = c(-1,36),
                             y.rev = TRUE)
addZone(tl09_stratplot, upper = c(27.5, 31.5, 32.5), lty = 1, lwd = 4, col = alpha("red", 0.2))
dev.off()

#######tl08 diatom supplement#########################

tl08_diatom_supp <- tl08_diatom
colnames(tl08_diatom_supp)<-gsub("[(M)]","",colnames(tl08_diatom_supp))
colnames(tl08_diatom_supp)<-gsub("B ","",colnames(tl08_diatom_supp))
colnames(tl08_diatom_supp)<-gsub("E ","",colnames(tl08_diatom_supp))
colnames(tl08_diatom_supp)<-gsub("F ","",colnames(tl08_diatom_supp))
colnames(tl08_diatom_supp)<-gsub("I ","",colnames(tl08_diatom_supp))

tl08_depth <- tl08_diatom_supp$depth
tl08_spec <- tl08_diatom_supp[,-1]
tl08_perc <- tl08_spec/rowSums(tl08_spec, na.rm = TRUE)*100
tl08_perc[is.na(tl08_perc)] <- 0

tl08_colsums <- colSums(tl08_perc)
tl08_sorted <- tl08_colsums[order(tl08_colsums, decreasing = TRUE)]
tl08_35_prep <- tl08_sorted[1:35]
tl08_nm <- names(tl08_35_prep)
tl08_pattern <- str_c(tl08_nm, collapse = "|")
tl08_35 <- select(tl08_perc, matches(tl08_pattern))

svg("figures/tl08_diatom_diagram.svg",
    width = 12.5,
    height = 5,
    bg = "white")
tl08_stratplot <- strat.plot(tl08_35,
                             cex.ylabel = 0.8,
                             title = "(D) Frenchmans Pond, TL08-2",
                             cex.title = 1,
                             ylabel = "Depth (cm)",
                             yvar = tl08_depth,
                             cex.xlabel = 0.8,
                             xRight = 0.95,
                             yTop = 0.6,
                             yBottom = 0.1,
                             srt.xlabel = 45,
                             y.tks = seq(0, 34, by = 2),
                             scale.percent = TRUE,
                             scale.minmax = TRUE,
                             plot.bar = TRUE,
                             plot.poly = TRUE,
                             plot.line = FALSE,
                             lwd.bar = 5,
                             col.bar = c(rep("darkred", 1),
                                         rep("darkblue", 4),
                                         rep("darkgray", 1),
                                         rep("darkgreen", 29),
                                         rep("black", 1)),
                             ylim = c(-1,34),
                             y.rev = TRUE)
addZone(tl08_stratplot, upper = 27.75, lty = 1, lwd = 4, col = alpha("red", 0.2))
dev.off()

###########tl18 PCA####
tl18_diatom_pca <- rda(sqrt(tl18_diatom_perc)) 

tl18_diatom_fort <- fortify(tl18_diatom_pca, axes = c(1,2), scaling = "sites")

tl18_diatom_sites <- tl18_diatom_fort[tl18_diatom_fort$Score %in% "sites",] %>% 
  mutate(depth = tl18_diatom$depth)
tl18_diatom_sp <- tl18_diatom_fort[tl18_diatom_fort$Score %in% "species",]

tl18_diatom_inertcomp <- inertcomp(tl18_diatom_pca) #contribution of each taxon to the total inertia
tl18_diatom_sel_sp <- tibble(taxon = rownames(tl18_diatom_inertcomp), inertcomp = tl18_diatom_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
tl18_diatom_sp_red <- tl18_diatom_sp[which(tl18_diatom_sp$Label %in% tl18_diatom_sel_sp$taxon), ] #leave 10 taxa with the highest contrubution to the total inertia

tl18_diatom_ve_prep <- tl18_diatom_pca$CA$eig / tl18_diatom_pca$tot.chi * 100
(tl18_diatom_PC1_ve <- round(((tl18_diatom_ve_prep / sum(tl18_diatom_ve_prep))[c(1)]) * 100, digits = 1))#25.9% expl. var.
(tl18_diatom_PC2_ve <- round(((tl18_diatom_ve_prep / sum(tl18_diatom_ve_prep))[c(2)]) * 100, digits = 1))#11.1% expl. var.

tl18_diatom_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", tl18_diatom_PC2_ve, "%)", sep = ""), x = paste("PC1 (", tl18_diatom_PC1_ve, "%)", sep = ""), title = "Borad Pond, TL18-2 diatom PCA plot") +
  geom_segment(data = tl18_diatom_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = tl18_diatom_sites, aes(x = PC1, y = PC2, color = depth)) +
  scale_color_viridis_c() +
  geom_point(data = filter(tl18_diatom_sites, depth %in% c(5.75, 23.75, 31.75, 35.75, 40.25)), aes(x = PC1, y = PC2), color = "red") +
  ggrepel::geom_text_repel(data = tl18_diatom_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(depth))) +
  ggrepel::geom_text_repel(data = tl18_diatom_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = tl18_diatom_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', size = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

ggsave(filename="figures/tl18_diatom_pca.svg", 
       plot = tl18_diatom_pca_plot, 
       device = svg, 
       width = 7, 
       height = 7, 
       units = "in")
