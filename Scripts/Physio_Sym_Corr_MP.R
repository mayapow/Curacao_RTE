####Physiology & symbiont ID correlations for Curaçao RTE 2026
#Maya Powell
#Last updated: Aug 7th 2026


####### Load libraries #########
library(lubridate)
library(tidyverse)
library(here)
library(car)
library(dplyr)
library(ggplot2)
library(powerjoin)
library(forcats)
library(emmeans)
library(see)
library(performance)
library(purrr)
library(rlang)
library(tibble)
library(psych)
library(corrplot)
library(ggpubr)
#library(vegan)
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#library(pairwiseAdonis)
library(lme4)
library(lmerTest) # note: this masks lme4::lmer() with its own version
#step() needs a model fit with lmerTest::lmer(), not lme4::lmer(), or it won't work

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
# write.csv(sym_physio, file = here("Data/sym_physio.csv"))
# manually filled in NAs in metadata to solve issues here before reading back in sym_physio

#read back in combined data
sym_physio <- read_csv(here("Data", "sym_physio.csv"))
sym_physio <- sym_physio %>% #and make it nice
  mutate(
    coral_species = factor(coral_species),
    native_site = factor(native_site, levels = c("Reef", "Bay")),
    transplant_site = factor(transplant_site, levels = c("Reef", "Bay")),
    sym_dom = factor(sym_dom),
    timepoint = factor(timepoint, levels = c("T0", "T4", "T12")),
    parent = factor(parent),
    native_transplant = factor(native_transplant, levels = c("Initial", "Native", "Transplant")),
    transplant_group = factor(
      transplant_group,
      levels = c("Reef Initial", "Bay Initial",
                 "Reef Native", "Reef Transplant",
                 "Bay Native", "Bay Transplant")))

##### Exploratory plots #####
#chlc2	chla	chlc2_sym	chla_sym	dry_weight	sym_density	calcification	pct_change_dw	pct_change_sa 
#create a function to plot all!!

#species
ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()
species_data <- list("Siderastrea siderea" = ss_data,"Porites sp." = pp_data)

#variables to examine
response_vars <- c("chla", "chla_sym", "dry_weight", "sym_density","calcification", "pct_change_dw", "pct_change_sa",
                   "PNET", "RDARK", "LEDR", "PGROSS", "PGROSS.RDAY")

response_labels <- list(
  chla = "Chlorophyll a (µg/cm2)",
  chla_sym = "Chlorophyll a per symbiont (pg/cell)",
  dry_weight = "Tissue biomass (mg/cm2)",
  sym_density = "Symbiont density (x10^6 cells/cm2)",
  calcification = "Calcification (mg/cm2/d)",
  pct_change_dw = "% change tissue biomass (T0-T12)",
  pct_change_sa = "% change surface area (T0-T12)",
  PNET = "Net Photosynthesis (µmol O2/cm2/hr)",
  RDARK = "Dark Respiration (µmol O2/cm2/hr)",
  LEDR = "Light Enhanced Dark Respiration (µmol O2/cm2/hr)",
  PGROSS = "Gross Photosynthesis (µmol O2/cm2/hr)",
  PGROSS.RDAY = "P:R Ratio"
)

#its2 type colors
sym_dom_colors = c("A4" = "#ffaabb","C1" = "#eedd88","C3" = "#77aadd", "C_other" = "#99ddff","D1" = "#994455")
site_colors = c("Bay Initial" = "coral4", "Bay Native" = "coral4", "Bay Transplant" = "salmon2", "Reef Transplant" = "lightblue3", "Reef Initial" = "dodgerblue4", "Reef Native" = "dodgerblue4")

#create a directory if needed to put output in
dir.create(here("Output", "Physiology"), recursive = TRUE, showWarnings = FALSE)
#dir.create(here("Output", "Tables"), recursive = TRUE, showWarnings = FALSE)

plot_physio_by_sym <- function(data, response) {
  ylab_txt <- response_labels[[response]] %||% response
  ggplot(data, aes(x = transplant_group, y = .data[[response]], color = sym_dom)) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
    scale_color_manual(values = sym_dom_colors) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(y = ylab_txt, color = "Dominant\nITS2 type") +
    facet_wrap(~ timepoint, scales = "free")
}

iwalk(species_data, function(d, sp) {
  walk(response_vars, function(resp) {
    p <- plot_physio_by_sym(d, resp)
    ggsave(
      here("Output", "Physiology", paste0(sp, "_", resp, "_sym.pdf")),
      p, height = 5, width = 12
    )
  })
})

#View each by changing variables: physio_plots[["dry_weight"]]


##### MODELS ######

#lmerTest::step() does backward elimination, so start from the most complex model

response_vars_with_time <- c("chla", "chla_sym", "dry_weight", "sym_density",
                             "PNET", "RDARK", "LEDR", "PGROSS", "PGROSS.RDAY")

predictor_names <- c("native_site", "transplant_site", "timepoint", "sym_dom")

clean_trait_data <- function(data, response_vars_with_time) {
  data %>%
    filter(
      !is.na(.data[[response]]), !is.na(sym_dom),
      !is.na(native_site), !is.na(transplant_site)
    ) %>%
    droplevels()
}

step_results_table <- tibble(response = character(), species = character(),
                             n_obs = integer(), final_formula = character())

for (response_name in response_vars_with_time) {
  for (species_name in names(species_data)) {
    
    this_data <- clean_trait_data(species_data[[species_name]], response_name, predictor_names)
    if (nrow(this_data) < 10) next
    
    model_formula <- as.formula(paste(
      response_name,
      "~ native_site * transplant_site * timepoint * sym_dom * Shannon * Observed * even + (1 | native_site/parent)"
    ))
    
    full_model <- tryCatch(lmerTest::lmer(model_formula, data = this_data, REML = TRUE),
                           error = function(e) NULL)
    if (is.null(full_model)) next
    
    step_result <- tryCatch(lmerTest::step(full_model), error = function(e) NULL)
    if (is.null(step_result)) next
    
    print(step_result)
    final_model <- get_model(step_result)
    print(summary(final_model))
    
    step_results_table <- step_results_table %>%
      add_row(response = response_name, species = species_name,
              n_obs = nrow(this_data),
              final_formula = paste(deparse(formula(final_model)), collapse = " "))
  }
}

print(step_results_table)
write_csv(step_results_table, here("Output", "Simple_Models", "step_selection_summary.csv"))
