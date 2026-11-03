####### RTE Physio + Symbiont Type Integration #####
#Maya Powell
#Last updated August 24th 2026

####### Load libraries #########
library(lubridate)
library(tidyverse)
library(here)
library(car)
library(dplyr)
library(ggplot2)
library(powerjoin)
library(forcats)
library(remotes)
library(emmeans)
#packageVersion("emmeans")
library(multcomp)
library(multcompView)
library(see)
library(performance)
library(purrr)
library(rlang)
#library(tibble)
#library(psych)
#library(corrplot)
library(ggpubr)
#library(vegan)
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#library(pairwiseAdonis)
library(lme4)
library(lmerTest) # note: this masks lme4::lmer() with its own version
#step() needs a model fit with lmerTest::lmer(), not lme4::lmer(), or it won't work
#emm_example("cld-multcomp") for some reason not working to do cld


###### Read in Data ########
# #sym data dominant type + metadata
# sym_dom <- read_csv(here("Data/ITS2.dominanttype.CW.RTE.csv"))
# #bacterial alpha diversity data
# micro_alpha <- read.csv(here("Data/CW.RTE.diversity.coral.ps.less.csv"))
# micro_alpha <- micro_alpha %>% select(id,Observed,Shannon,even)
# sym_micro <- sym_dom %>% left_join(micro_alpha, by = "id")
# 
# #physiology from sarah
# all_physio <- read_csv(here("Data/all_physio_notransplant.csv"))
# all_physio <- all_physio %>% select(sample_id, chlc2,	chla,	chlc2_sym,	chla_sym,	dry_weight,	sym_density)
# calc_data <- read_csv(here("Data/calc_data.csv"))
# calc_data <- calc_data %>% select(sample_id, calcification,	pct_change_dw,	pct_change_sa)
# 
# #join physio & sym data
# all_physio <- all_physio %>% left_join(calc_data, by = "sample_id")
# sym_physio <- all_physio %>% left_join(sym_micro, by = "sample_id")
# 
# metabolic_data <- read.csv(here("Data/metabolic_data.csv"))
# metabolic_data <- metabolic_data %>% select(sample_id,PNET,RDARK,LEDR,PGROSS,PGROSS.RDAY)
# sym_physio <- sym_physio %>% left_join(metabolic_data, by = "sample_id")
# sym_physio$...1 <- NULL #remove first row bc it's annoying
# write.csv(sym_physio, file = here("Data/sym_physio.csv"), row.names=FALSE)
# manually filled in NAs in metadata to solve issues here before reading back in sym_physio

#read back in combined data
sym_physio <- read_csv(here("Data", "sym_physio.csv"), locale=locale(encoding="latin1")) #had to add this encoding to be able to see in r for some reason something weird about a character width
sym_physio <- sym_physio %>% #and make it nice/full of factors
  mutate(
    coral_species = factor(coral_species),
    native_site = factor(native_site, levels = c("Reef", "Bay")),
    native_transplant = factor(native_transplant, levels = c("Native", "Transplant")),
    sym_dom = factor(sym_dom),
    timepoint = factor(timepoint, levels = c("T0", "T4", "T12")),
    parent = factor(parent),
    native_transplant = factor(native_transplant, levels = c("Initial", "Native", "Transplant")),
    transplant_group_sampled = factor(transplant_group_sampled,
                                      levels = c("Bay to Bay", "Bay to Reef", "Reef to Bay", "Reef to Reef")),
    transplant_group = factor(
      transplant_group,
      levels = c("Reef Initial", "Bay Initial","Reef Native", "Reef Transplant", "Bay Native", "Bay Transplant")))

#complete data transformations for all variables based on Sarah data (see info from Sarah models at bottom of script)
#this will create normal data for all variables to use in models (except few that use glmer transformations)
sym_physio <- sym_physio %>%
  mutate(
    chla_log = log(chla),
    chla_sqrt = sqrt(chla),
    chla_sym_log = log(chla_sym),
    chla_sym_sqrt = sqrt(chla_sym),
    sym_density_log = log(sym_density),
    sym_density_sqrt = sqrt(sym_density),
    PGROSS_sqrt = sqrt(PGROSS),
    RDARK_abs = abs(RDARK),
    LEDR_abs = abs(LEDR),
    LEDR_abs_log = log(LEDR_abs)
  )
#manually removed outliers from sym_physio for variables based on Sarah's findings (see at bottom of script, only 3 from calcification)


##### Simple plots of all variables to look at all data #####

##### Exploratory plots #####
#chlc2	chla	chlc2_sym	chla_sym	dry_weight	sym_density	calcification	pct_change_dw	pct_change_sa 
#create a function to plot all!!

#species
ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()
species_data <- list("Siderastrea siderea" = ss_data,"Porites sp." = pp_data)

#variables to examine
response_vars_plots <- c("chla", "chla_sym", "dry_weight", "sym_density","calcification",
                         "PNET", "RDARK", "LEDR", "PGROSS", "PGROSS.RDAY")

response_labels <- list(
  chla = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})),
  chla_sym = expression("Chlorophyll a" ~ (pg ~ cell^{-1})),
  dry_weight = expression("Biomass" ~ (mg ~ cm^{-2})),
  sym_density = expression("Sym" ~ (x10^{6} ~ cells ~ cm^{-2})),
  calcification = expression("Calcification" ~ (mg ~ cm^{-2} ~ d^{-1})),
  PNET = expression("NP" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  RDARK = expression("R" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  LEDR = expression("LEDR" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  PGROSS = expression("GP" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  PGROSS.RDAY = "P:R Ratio")

#its2 type colors
sym_dom_colors = c("A4" = "#ffaabb","C1" = "#eedd88","C3" = "#77aadd", "C_other" = "#99ddff","D1" = "#994455")
site_colors = c("Bay Initial" = "coral4", "Bay Native" = "coral4", "Bay Transplant" = "salmon2", "Reef Transplant" = "lightblue3", "Reef Initial" = "dodgerblue4", "Reef Native" = "dodgerblue4")

#create a directory if needed to put output in
dir.create(here("Output", "Physiology", "Plots"), recursive = TRUE, showWarnings = FALSE)

#standard error calculation
se_fun <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

#get means grouped by site
physio_mean_se <- function(data, response) {
  data %>%
    group_by(transplant_group, timepoint) %>%
    summarise(
      n    = sum(!is.na(.data[[response]])),
      mean = mean(.data[[response]], na.rm = TRUE),
      se   = se_fun(.data[[response]]),
      .groups = "drop"
    )
}

###Physiology Plots by Sym Time Site 
#create plots with means as well as jittered points
plot_physio_site_time_sym <- function(data, response) {
  ylab_txt <- response_labels[[response]] %||% response
  summary_data <- physio_mean_se(data, response)
  
  ggplot(data, aes(x = transplant_group, y = .data[[response]])) +
    geom_jitter(aes(color = sym_dom), width = 0.15, alpha = 0.5, size = 2) +
    geom_errorbar(data = summary_data,
                  aes(x = transplant_group, ymin = mean - se, ymax = mean + se),
                  width = 0.2, linewidth = 0.7, inherit.aes = FALSE) +
    geom_point(data = summary_data,
               aes(x = transplant_group, y = mean, fill = transplant_group),
               shape = 21, color = "black", size = 4, inherit.aes = FALSE) +
    scale_color_manual(values = sym_dom_colors) +
    scale_fill_manual(values = site_colors, guide = "none") +
    theme_bw(base_size = 22) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(y = ylab_txt, color = "Dominant\nITS2 type") +
    facet_wrap(~ timepoint, scales = "free_x")
}

#and run plots!

# #Physio site time sym
# iwalk(species_data, function(d, sp) {
#   walk(response_vars_plots, function(resp) {
#     p <- plot_physio_site_time_sym(d, resp)
#     ggsave(
#       here("Output", "Physiology", "Plots", paste0(sp, "_", resp, "site_time_sym.pdf")),
#       p, height = 5, width = 12
#     )
#   })
# })

##### Models of symbiont type and pairwise comparisons to add to plots #####

###Siderastrea siderea
#"chla_log", "chla_sym_log", "sym_density_log", "dry_weight", 
#"PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs_log", "calcification"
sym_physio_sid <- sym_physio %>% filter(coral_species == "Siderastrea siderea")
sym_physio_por <- sym_physio %>% filter(coral_species == "Porites sp.")

###all models across time: models with a significant effect of symbiont identity to examine further
##only affects bulk physiology
sid_mod_PR <- lmer(PGROSS.RDAY	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_PR)
sid_mod_PNET <- lmer(PNET	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_PNET)
sid_mod_PGROSS <- lmer(PGROSS	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_PGROSS)
sid_mod_RDARK <- lmer(RDARK_abs	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_RDARK)
sid_mod_LEDR <- lmer(LEDR_abs_log	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_LEDR)
sid_mod_chla <- lmer(chla_log	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_chla)
sid_mod_chla_sym <- lmer(chla_sym_log	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_chla_sym)
sid_mod_sym <- lmer(sym_density_log	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym)
sid_mod_dw <- lmer(dry_weight	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_dw)
sid_mod_calc <- lmer(calcification	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_calc)

summary(sid_mod_chla)
ranova(sid_mod_chla) #ns

sid_emm_chla <-emmeans(sid_mod_chla, pairwise ~  native_site * native_transplant) 
sid_emm_chla <- emmeans(sid_mod_chla, specs = c("sym_dom", "native_site"))
pairs(sid_emm_chla)

###Models only including symbiont identity
#PR
sid_mod_sym_PR <- lmer(PGROSS.RDAY	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_PR) #sig
sid_emm_PR <- emmeans::emmeans(sid_mod_sym_PR, ~ sym_dom)
pairs(sid_emm_PR)
#PNET
sid_mod_sym_PNET <- lmer(PNET	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_PNET) #sig
sid_emm_PNET <- emmeans::emmeans(sid_mod_sym_PNET, ~ sym_dom)
pairs(sid_emm_PNET)
#PGROSS
sid_mod_sym_PGROSS <- lmer(PGROSS	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_PGROSS) #sig
sid_emm_PGROSS <- emmeans::emmeans(sid_mod_sym_PGROSS, ~ sym_dom)
pairs(sid_emm_PGROSS)
#RDARK abs
sid_mod_sym_RDARK <- lmer(RDARK_abs	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_RDARK) #ns
#LEDR abs log
sid_mod_sym_LEDR <- lmer(LEDR_abs_log	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_LEDR) #ns
#chla log
sid_mod_sym_chla <- lmer(chla_log	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_chla) #sig
sid_emm_chla <- emmeans::emmeans(sid_mod_sym_chla, ~ sym_dom)
pairs(sid_emm_chla)
#chla sym log
sid_mod_sym_chla_sym <- lmer(chla_sym_log	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_chla_sym) #sig
sid_emm_chla_sym <- emmeans::emmeans(sid_mod_sym_chla_sym, ~ sym_dom)
pairs(sid_emm_chla_sym)
#sym density
sid_mod_sym_sym <- lmer(sym_density_log	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_sym) #sig
sid_emm_sym <- emmeans::emmeans(sid_mod_sym_sym, ~ sym_dom)
pairs(sid_emm_sym)
#dry weight
sid_mod_sym_dw <- lmer(dry_weight	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_dw) #ns
#calcification
sid_mod_sym_calc <- lmer(calcification	~  sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(sid_mod_sym_calc) #sig
sid_emm_calc <- emmeans::emmeans(sid_mod_sym_calc, ~ sym_dom)
pairs(sid_emm_calc)

### Branching Porites sp.
#"chla_sqrt", "chla_sym_sqrt", "sym_density_sqrt", "dry_weight", "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs", "calcification"
#normality of all models checked using the check_model from above which prints for best model
#can be found in Output, Physiology

###Models only including symbiont identity
#PR
por_mod_sym_PR <- lmer(PGROSS.RDAY	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_PR) #ns
#PNET
por_mod_sym_PNET <- lmer(PNET	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_PNET) #ns
#PGROSS
por_mod_sym_PGROSS <- lmer(PGROSS	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_PGROSS) #ns
#RDARK_abs
por_mod_sym_RDARK <- lmer(RDARK_abs	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_RDARK) #ns
#LEDR_abs
por_mod_sym_LEDR <- lmer(LEDR_abs	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_LEDR) #ns
#chla sqrt
por_mod_sym_chla <- lmer(chla_sqrt	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_chla) #ns
#chla sym sqrt
por_mod_sym_chla_sym <- lmer(chla_sym_sqrt	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_chla_sym) #sig
por_emm_chla_sym <- emmeans::emmeans(por_mod_sym_chla_sym, ~ sym_dom)
pairs(por_emm_chla_sym)
#sym density sqrt
por_mod_sym_sym <- lmer(sym_density_sqrt	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_sym) #sig
por_emm_sym <- emmeans::emmeans(por_mod_sym_sym, ~ sym_dom)
pairs(por_emm_sym)
#dry weight
por_mod_sym_dw <- lmer(dry_weight	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_dw) #sig
por_emm_dw <- emmeans::emmeans(por_mod_sym_dw, ~ sym_dom)
pairs(por_emm_dw)
#calcification
por_mod_sym_calc <- lmer(calcification	~  sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
anova(por_mod_sym_calc) #sig
por_emm_calc <- emmeans::emmeans(por_mod_sym_calc, ~ sym_dom)
pairs(por_emm_calc)

###Physiology Plots by Sym
sym_mean_se <- function(data, response) {
  data %>%
    group_by(sym_dom) %>%
    summarise(
      n    = sum(!is.na(.data[[response]])),
      mean = mean(.data[[response]], na.rm = TRUE),
      se   = se_fun(.data[[response]]),
      .groups = "drop"
    )
}

#create plots with means as well as jittered points
plot_physio_sym <- function(data, response) {
  ylab_txt <- response_labels[[response]] %||% response
  summary_data <- sym_mean_se(data, response)
  
  ggplot(data, aes(x = sym_dom, y = .data[[response]])) +
    geom_jitter(aes(color = sym_dom), width = 0.15, alpha = 0.8, size = 2) +
    geom_errorbar(data = summary_data,
                  aes(x = sym_dom, ymin = mean - se, ymax = mean + se),
                  width = 0.2, linewidth = 0.7, inherit.aes = FALSE) +
    geom_point(data = summary_data,
               aes(x = sym_dom, y = mean),
               fill = "black", size = 4, inherit.aes = FALSE) +
    scale_color_manual(values = sym_dom_colors) +
    theme_bw(base_size = 22) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(y = ylab_txt, color = "Dominant\nITS2 type")
}

#Physio sym
physio_sym_plots <- imap(species_data, function(d, sp) {
  set_names(
    map(response_vars_plots, function(resp) plot_physio_sym(d, resp)),
    response_vars_plots
  )
})

# #Save all plots
# iwalk(physio_sym_plots, function(species_plots, sp) {
#   iwalk(species_plots, function(p, resp) {
#     ggsave(
#       here("Output", "Physiology", "Plots", paste0(sp, "_", resp, "_sym.pdf")),
#       p, height = 5, width = 8
#     )
#   })
# })

#examine each plot:
#physio_sym_plots[["Siderastrea siderea"]][["chla"]]

#Update and annotate plots with letters for significance
#"chla", "chla_sym", "dry_weight", "sym_density","calcification",
#"PNET", "RDARK", "LEDR", "PGROSS", "PGROSS.RDAY"

#siderastrea siderea
sid_chla_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["chla"]] +
  ylim(3,25) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("ab", "a", "b"))
sid_chla_sym_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["chla_sym"]] +
  ylim(2,14) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "b"))
sid_dw_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["dry_weight"]]
sid_sym_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["sym_density"]] +
  ylim(0,5) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "b"))
sid_calc_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["calcification"]] +
  ylim(0,2.5) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "b"))
sid_PNET_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["PNET"]] +
  ylim(0,3) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "ab"))
sid_PGROSS_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["PGROSS"]] +
  ylim(0.5,4) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "ab"))
sid_RDARK_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["RDARK"]]
sid_LEDR_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["LEDR"]]
sid_PR_sym_plot <- physio_sym_plots[["Siderastrea siderea"]][["PGROSS.RDAY"]] +
  ylim(0.25,2) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "b"))

sid_all_physio_sym <- ggarrange(
  sid_chla_sym_plot,sid_chla_sym_sym_plot,sid_sym_sym_plot,sid_dw_sym_plot,sid_calc_sym_plot,
  sid_PR_sym_plot,sid_PNET_sym_plot,sid_PGROSS_sym_plot,sid_RDARK_sym_plot,sid_LEDR_sym_plot,
  ncol = 5, nrow = 2, legend = "right", common.legend = TRUE, 
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I","J"),
  font.label = list(size = 30, color = "black", face = "bold", family = NULL)
  )

ggsave(here("Output/Physiology/Plots/sid_all_physio_sym.pdf"),
       sid_all_physio_sym, height = 10, width = 20)

#branching porites sp
por_chla_sym_plot <- physio_sym_plots[["Porites sp."]][["chla"]]
por_chla_sym_sym_plot <- physio_sym_plots[["Porites sp."]][["chla_sym"]] +
  ylim(2,17) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "ab"))
por_dw_sym_plot <- physio_sym_plots[["Porites sp."]][["dry_weight"]] +
  ylim(0,42) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "a"))
por_sym_sym_plot <- physio_sym_plots[["Porites sp."]][["sym_density"]]
por_calc_sym_plot <- physio_sym_plots[["Porites sp."]][["calcification"]] +
  ylim(0.5,3) +
  stat_summary(geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "ab"))
por_PNET_sym_plot <- physio_sym_plots[["Porites sp."]][["PNET"]]
por_PGROSS_sym_plot <- physio_sym_plots[["Porites sp."]][["PGROSS"]]
por_RDARK_sym_plot <- physio_sym_plots[["Porites sp."]][["RDARK"]]
por_LEDR_sym_plot <- physio_sym_plots[["Porites sp."]][["LEDR"]]
por_PR_sym_plot <- physio_sym_plots[["Porites sp."]][["PGROSS.RDAY"]]

por_all_physio_sym <- ggarrange(
  por_chla_sym_plot,por_chla_sym_sym_plot,por_sym_sym_plot,por_dw_sym_plot,por_calc_sym_plot,
  por_PR_sym_plot,por_PNET_sym_plot,por_PGROSS_sym_plot,por_RDARK_sym_plot,por_LEDR_sym_plot,
  ncol = 5, nrow = 2, legend = "right", common.legend = TRUE, 
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I","J"),
  font.label = list(size = 30, color = "black", face = "bold", family = NULL)
)

ggsave(here("Output/Physiology/Plots/por_all_physio_sym.pdf"),
       por_all_physio_sym, height = 10, width = 20)


####### Models and info used by Sarah for each variable for reference #######

# #chlorophyll a per surface area
# #no outliers removed
# #sid log transformed data for normality
# model1_SID_Tchla<-lmer(SID_log_Tchla ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = SID) #use this one 
# #por sqrt transformed for normality
# model1_POR_Tchla<-lmer(POR_sqrt_Tchla ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = POR) #use this one 
# 
# #chlorophyll a per cell
# #no outliers removed
# #sid log transformed for normality
# model1_SID_TchlaCell<-lmer(SID_log_TchlaCell ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = SID) #use this one 
# #por sqrt transformed for normality
# #with random effect - singular fit can take out
# model1_POR_TchlaCell<-lm(POR_sqrt_Tchla ~ Timepoint*Treatment*Origin, + (1| Origin:`Unique Genotype`), data = POR) #use this one 
# 
# #symbiont density
# #no outliers removed
# #sid log transformed for normality
# model1_SID_ZOOXS<-lmer(SID_log_ZOOXS ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = SID) #use this one 
# #por sqrt transformed for normality
# model1_POR_ZOOXS<-lmer(POR_sqrt_ZOOXS ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = POR) #use this one
# 
# #biomass
# #no outliers removed
# #no transformations
# model1_SID_Biomass<-lmer(DTW ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = SID)
# model1_POR_Biomass<-lmer(DTW ~ Timepoint*Treatment*Origin + (1| Origin:`Unique Genotype`), data = POR)
# 
# #biomass at T12 for reaction norms
# #no transformations
# #por no outliers removed
# #sid remove SID SPNB3 at T12
# SID_BiomassT12_outliers<-SID_biomass_T12 %>%
#   identify_outliers(DTW) #one extreme outlier SID SPNB 3
# model1_SID_Biomass<-lmer(DTW ~ Treatment*Origin + (1| Origin:`Unique Genotype`), data = SID_biomass_T12) 
# model1_POR_Biomass<-lmer(DTW ~ Treatment*Origin + (1| Origin:`Unique Genotype`), data = POR_biomass_T12)
# 
# #P:R ratio (PGROSSRDAY)
# #no outliers removed
# #no transformations
# model1_Sid_PGROSSRDAY<-lmer(PGROSSRDAY ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = Incubation_SID)
# model1_Por_PGROSSRDAY<-glmer(PGROSSRDAY ~ Timepoint*Treatment*Origin + (1|Origin:Genotype), data = Incubation_POR, family = Gamma(link = "identity"))
# 
# #Net photosynthesis (PNET)
# #no outliers removed
# #no transformations
# model1_Sid_PNET<-lmer(PNET ~ Timepoint*Treatment*Origin + (1|Origin:Genotype), data = Incubation_SID)
# model1_Por_PNET<-lmer(PNET ~ Timepoint*Treatment*Origin + (1|Origin:Genotype), data = Incubation_POR)
# 
# #Gross photosynthesis (PGROSS)
# #no outliers
# #sid sqrt transformed for normality
# model1_Sid_PGROSS<-lmer(SID_sqrt_PGROSS ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = Incubation_SID) #use this one
# model1_Por_PGROSS<-lmer(PGROSS ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = Incubation_POR)
# 
# #Respiration (RDARK)
# #no outliers removed
# #no transformations
# #sarah uses model without random effect for POR
# #for sid: specified an optimizer called "bobyqa". Different optimizers result in same or similar Anova results
# model1_Sid_RDARK<-glmer(SID_abs_RDARK~ Timepoint*Treatment*Origin + (1|Origin:Genotype), family=Gamma(link= "log"), control = glmerControl(optimizer = "bobyqa"), data = Incubation_SID)
# model1_Por_RDARK<-glmer(abs(RDARK) ~ Timepoint*Treatment*Origin + (1|Origin:Genotype), data = Incubation_POR, family = Gamma(link = "log"))
# 
# #Light enhanced dark respiration (LEDR)
# #no outliers removed
# #sid log transformed for normality
# model1_Sid_LEDR<-lmer(SID_log_LEDR ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = Incubation_SID) #use this one
# model1_Por_LEDR<-lmer(LEDR ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = Incubation_POR)
# 
# #Calcification
# #no transformations
# #sid outliers to remove in BWdata: BWdata_SID2<-BWdata_SID[-c(51),] #ID SID SPNA2 (from October 2022)
# model1.1_Sid<-lmer(Calcification ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = BWdata_SID2)#use this one
# #por outliers to remove in BWdata: BWdata_POR2<-BWdata_POR[-c(55,59),]
# model2_Por<-lmer(Calcification ~ Timepoint*Treatment*Origin + (1| Origin:Genotype), data = BWdata_POR2)#use this one
# 
# #Calcification at T12
# #same outliers removed as in calcification por and sid
# model2_Por_T12<-lmer(Calcification ~ Treatment*Origin + (1| Origin:Genotype), data = BWdata_POR2_T12)#use this one
# model2_SID_T12<-lmer(Calcification ~ Treatment*Origin + (1| Origin:Genotype), data = BWdata_SID2_T12)
