####### Chlorophyll a model comparison - Curaçao RTE 2026 #######
# Maya Powell
# Adapted with Claude: Aug 19 2026
#
# Simplified, single-trait version of the larger symbiont/physiology
# script, focused only on chlorophyll a (chla). Differences from the
# big script:
#   - Only chla is modeled (no trait loop).
#   - Every candidate model is written out explicitly as its own named
#     lmer() call -- no boolean/tribble-driven formula generation.
#   - Random effect is `parent` NESTED WITHIN `native_site`
#     (1 | native_site/parent) instead of a plain (1|parent).
#   - Winning model per species is checked with performance::check_model()
#     before being trusted.
#   - Final plot is a plain jittered scatter (chla over timepoint, faceted
#     by transplant group, colored by dominant symbiont type) -- no
#     boxplots, no compact letter displays.
#
# Section 0: libraries, colors, data read-in, species split
# Section 1: explicit chla model set (one function per species' data)
# Section 2: fit all models, per species
# Section 3: AIC comparison table
# Section 4: check_model() diagnostics on each species' winning model
# Section 5: refit winner with REML, Anova + R2
# Section 6: simple plot - chla x timepoint, facet by site, color by sym_dom


####### Load libraries #########
library(tidyverse)
library(here)
library(car)
library(lme4)
library(emmeans)
library(see)
library(performance)   # check_model(), r2() -- check_model()'s plot output
                        # depends on the 'see' package; install.packages("see")
                        # if the diagnostic PDF comes out blank

####### Color palette #########
sym_dom_colors <- c("A4" = "#ffaabb", "C1" = "#eedd88", "C3" = "#77aadd",
                     "C_other" = "#99ddff", "D1" = "#994455")

###### Read in data ########
sym_physio <- read_csv(here("Data", "sym_physio.csv"),
                        locale = locale(encoding = "latin1"))

###### Clean & prep data ########

sym_physio <- sym_physio %>%
  mutate(
    coral_species   = factor(coral_species),
    native_site     = factor(native_site, levels = c("Reef", "Bay")),
    transplant_site = factor(transplant_site, levels = c("Reef", "Bay")),
    sym_dom         = factor(sym_dom),
    timepoint       = factor(timepoint, levels = c("T0", "T4", "T12")),
    parent          = factor(parent),
    transplant_group = factor(
      transplant_group,
      levels = c("Reef Initial", "Bay Initial",
                 "Reef Native", "Reef Transplant",
                 "Bay Native", "Bay Transplant")
    )
  ) %>%
  filter(!is.na(chla), !is.na(sym_dom), !is.na(native_site), !is.na(transplant_site)) %>%
  droplevels()

ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()

dir.create(here("Output", "Chla_Model_Comparison"), recursive = TRUE, showWarnings = FALSE)

####### 1. Explicit chla model set ########
# Same conceptual candidate space as the big script (site structure x
# symbiont x time), but every formula is written out by hand so it's
# immediately readable/editable. Random effect throughout: parent nested
# within native_site, i.e. (1 | native_site/parent).
#
# fit_one() just wraps lmer() in a tryCatch so one non-converging model
# (most likely among the busiest sym_dom x site x time models near the
# bottom of the list, given how few samples some rarer symbiont types
# have) doesn't stop the whole script -- it will simply come back NULL
# and get dropped from the comparison table.

fit_one <- function(formula, data) {
  tryCatch(
    suppressWarnings(suppressMessages(
      lmer(formula, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    )),
    error = function(e) NULL
  )
}

fit_chla_models <- function(d) {
  list(

    null = fit_one(
      chla ~ timepoint + (1 | native_site/parent), d),

    # ---- site only ----
    native = fit_one(
      chla ~ native_site + timepoint + (1 | native_site/parent), d),
    native_tp = fit_one(
      chla ~ native_site * timepoint + (1 | native_site/parent), d),

    transplant = fit_one(
      chla ~ transplant_site + timepoint + (1 | native_site/parent), d),
    transplant_tp = fit_one(
      chla ~ transplant_site * timepoint + (1 | native_site/parent), d),

    site_add = fit_one(
      chla ~ native_site + transplant_site + timepoint + (1 | native_site/parent), d),
    site_add_tp = fit_one(
      chla ~ (native_site + transplant_site) * timepoint + (1 | native_site/parent), d),

    site_int = fit_one(
      chla ~ native_site * transplant_site + timepoint + (1 | native_site/parent), d),
    site_int_tp = fit_one(
      chla ~ native_site * transplant_site * timepoint + (1 | native_site/parent), d),

    # ---- symbiont only ----
    sym = fit_one(
      chla ~ sym_dom + timepoint + (1 | native_site/parent), d),
    sym_tp = fit_one(
      chla ~ sym_dom * timepoint + (1 | native_site/parent), d),

    # ---- site + symbiont, additive with each other ----
    native_sym = fit_one(
      chla ~ native_site + sym_dom + timepoint + (1 | native_site/parent), d),
    native_sym_tp = fit_one(
      chla ~ (native_site + sym_dom) * timepoint + (1 | native_site/parent), d),

    transplant_sym = fit_one(
      chla ~ transplant_site + sym_dom + timepoint + (1 | native_site/parent), d),
    transplant_sym_tp = fit_one(
      chla ~ (transplant_site + sym_dom) * timepoint + (1 | native_site/parent), d),

    site_add_sym_add = fit_one(
      chla ~ native_site + transplant_site + sym_dom + timepoint + (1 | native_site/parent), d),
    site_add_sym_add_tp = fit_one(
      chla ~ (native_site + transplant_site + sym_dom) * timepoint + (1 | native_site/parent), d),

    site_int_sym_add = fit_one(
      chla ~ native_site * transplant_site + sym_dom + timepoint + (1 | native_site/parent), d),
    site_int_sym_add_tp = fit_one(
      chla ~ (native_site * transplant_site + sym_dom) * timepoint + (1 | native_site/parent), d),

    # ---- site x symbiont, interactive with each other ----
    native_sym_int = fit_one(
      chla ~ native_site * sym_dom + timepoint + (1 | native_site/parent), d),
    native_sym_int_tp = fit_one(
      chla ~ native_site * sym_dom * timepoint + (1 | native_site/parent), d),

    transplant_sym_int = fit_one(
      chla ~ transplant_site * sym_dom + timepoint + (1 | native_site/parent), d),
    transplant_sym_int_tp = fit_one(
      chla ~ transplant_site * sym_dom * timepoint + (1 | native_site/parent), d),

    site_add_sym_int = fit_one(
      chla ~ (native_site + transplant_site) * sym_dom + timepoint + (1 | native_site/parent), d),
    site_add_sym_int_tp = fit_one(
      chla ~ (native_site + transplant_site) * sym_dom * timepoint + (1 | native_site/parent), d),

    site_int_sym_int = fit_one(
      chla ~ native_site * transplant_site * sym_dom + timepoint + (1 | native_site/parent), d),
    site_int_sym_int_tp = fit_one(
      chla ~ native_site * transplant_site * sym_dom * timepoint + (1 | native_site/parent), d)
  )
}

####### 2. Fit all models, per species ########

ss_models <- fit_chla_models(ss_data)
pp_models <- fit_chla_models(pp_data)

####### 3. AIC comparison table ########

aic_compare <- function(models_list, species) {
  models_list <- models_list[!map_lgl(models_list, is.null)]
  tibble(
    species  = species,
    model    = names(models_list),
    df       = map_dbl(models_list, ~ attr(logLik(.x), "df")),
    AIC      = map_dbl(models_list, AIC),
    singular = map_lgl(models_list, isSingular)
  ) %>%
    mutate(
      deltaAIC = AIC - min(AIC),
      weight   = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
    ) %>%
    arrange(AIC)
}

ss_aic <- aic_compare(ss_models, "Siderastrea siderea")
pp_aic <- aic_compare(pp_models, "Porites sp.")

chla_model_comparison <- bind_rows(ss_aic, pp_aic)
print(chla_model_comparison)

write_csv(chla_model_comparison, here("Output", "Chla_Model_Comparison", "chla_model_comparison.csv"))

best_model_name_ss <- ss_aic$model[1]
best_model_name_pp <- pp_aic$model[1]

cat("\nBest model for S. siderea chla:", best_model_name_ss, "\n")
cat("Best model for Porites chla:", best_model_name_pp, "\n")

####### 4. check_model() diagnostics on the winning model per species ########
# Diagnostics are run on the REML fit (standard practice for inference/
# diagnostics, vs. the ML fits used above for fair AIC comparison).
# Refit explicitly from the winning formula (rather than update()) so
# there's no dependence on the fit_chla_models() function's internal
# environment still being around.

refit_reml <- function(model, data) {
  lmer(formula(model), data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
}

best_mod_ss <- refit_reml(ss_models[[best_model_name_ss]], ss_data)
best_mod_pp <- refit_reml(pp_models[[best_model_name_pp]], pp_data)

pdf(here("Output", "Chla_Model_Comparison", "chla_checkmodel_Ssiderea.pdf"), width = 12, height = 10)
print(check_model(best_mod_ss))
dev.off()

pdf(here("Output", "Chla_Model_Comparison", "chla_checkmodel_Porites.pdf"), width = 12, height = 10)
print(check_model(best_mod_pp))
dev.off()

# View interactively instead of/in addition to saving, e.g.:
# check_model(best_mod_ss)

####### 5. Refit winner (REML) - Anova + R2 ########

cat("\n========== S. siderea: chla, best model =", best_model_name_ss, "==========\n")
print(car::Anova(best_mod_ss))
print(performance::r2(best_mod_ss))

cat("\n========== Porites sp.: chla, best model =", best_model_name_pp, "==========\n")
print(car::Anova(best_mod_pp))
print(performance::r2(best_mod_pp))

####### 6. Simple plot: chla over time, faceted by site, colored by symbiont ########
# Plain jittered points only -- no boxplots, no significance letters.

plot_chla_simple <- function(data, species_name) {
  ggplot(data, aes(x = timepoint, y = chla, color = sym_dom)) +
    geom_jitter(width = 0.15, alpha = 0.85, size = 2.5) +
    scale_color_manual(values = sym_dom_colors, name = "Dominant\nITS2 type") +
    facet_wrap(~ transplant_group, nrow = 1) +
    theme_bw(base_size = 15) +
    theme(axis.title.x = element_blank()) +
    labs(y = "Chlorophyll a (µg/cm2)", title = paste("Chlorophyll a -", species_name))
}

p_chla_ss <- plot_chla_simple(ss_data, "Siderastrea siderea")
p_chla_pp <- plot_chla_simple(pp_data, "Porites sp.")

ggsave(here("Output", "Chla_Model_Comparison", "chla_over_time_Ssiderea.pdf"),
       p_chla_ss, width = 12, height = 5)
ggsave(here("Output", "Chla_Model_Comparison", "chla_over_time_Porites.pdf"),
       p_chla_pp, width = 12, height = 5)

# View interactively, e.g.:
# p_chla_ss
# p_chla_pp
