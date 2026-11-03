####### Physiology model comparison - simple version for teaching #######
# Maya Powell
# Adapted with Claude: Aug 19 2026
#
# Restructured into clearly separate steps so individual pieces can be
# inspected/re-run on their own -- this is especially useful right now
# for figuring out why PNET/RDARK/LEDR/PGROSS/PGROSS.RDAY aren't fitting.
#
#  PART 1: PLOTS ONLY. Make and save the raw-data jitter plot for every
#  trait, for both species, before any modeling happens. Look at these
#  first -- if a trait's data looks obviously wrong (all one value, huge
#  outliers, all NA for one species, etc.) that alone may explain why
#  models won't fit.
#
#  PART 2: FIT + CHECK + AIC COMPARE. For every trait x species, fit all
#  candidate models (43 for the T0/T4/T12 traits -- the full additive /
#  partially interactive / fully interactive set across native_site,
#  transplant_site, sym_dom, and timepoint; 14 for the T12-only traits,
#  which have no timepoint at all), compare them by AIC, and run check_model() on the
#  best one. Every fitted model (or NULL, if it failed) is stored in
#  all_fitted_models, and every AIC table is stored in all_aic_tables, so
#  you can pull any individual piece back out afterward, e.g.:
#     all_fitted_models[["PNET"]][["Porites sp."]]$site_int
#     all_fitted_models[["PNET"]][["Porites sp."]]$site_int_time
#     all_aic_tables[["PNET"]][["Porites sp."]]
#
#  PART 3: WHICH PREDICTORS WIN. A single summary table, one row per
#  trait x species, built from the AIC tables saved in Part 2 -- no
#  refitting here, just pulling together what already happened.
#
# This version also changed fit_one(): instead of silently swallowing
# every warning and error, it now PRINTS them. That's more cluttered, but
# it should make it obvious exactly what's going wrong with the PNET-type
# variables.

####### Load libraries #########
library(tidyverse)
library(here)
library(car)
library(lme4)
library(see)
library(performance)   # check_model(), r2() -- check_model()'s plot output
                        # depends on the 'see' package; install.packages("see")
                        # if the diagnostic PDF comes out blank

####### Color palette #########
sym_dom_colors <- c("A4" = "#ffaabb", "C1" = "#eedd88", "C3" = "#77aadd",
                     "C_other" = "#99ddff", "D1" = "#994455")

###### Read in data ########
sym_physio <- read_csv(here("Data", "sym_physio.csv"))

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
  )

# Split into the two species.
ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()

species_data <- list(
  "Siderastrea siderea" = ss_data,
  "Porites sp."          = pp_data
)

# Traits measured at T0, T4, AND T12
response_vars_with_time <- c("chla", "chla_sym", "dry_weight", "sym_density",
                              "PNET", "RDARK", "LEDR", "PGROSS", "PGROSS.RDAY")

# Traits that are already "change from T0 to T12" -- T12 rows only
response_vars_no_time <- c("calcification", "pct_change_dw", "pct_change_sa")

# Make the output folder
dir.create(here("Output", "Simple_Models"), recursive = TRUE, showWarnings = FALSE)


####### Helper: safely fit one model, but PRINT what goes wrong #########

fit_one <- function(model_formula, data) {
  tryCatch({
    withCallingHandlers({
      lmer(model_formula, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    }, warning = function(w) {
      cat("    [warning]", deparse(model_formula), "->", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    cat("    [ERROR]  ", deparse(model_formula), "->", conditionMessage(e), "\n")
    NULL
  })
}


####### Build all model possibilities #########

fit_models_with_time <- function(response_name, data) {
  # ---- 0 predictors ----
  f_null                 <- as.formula(paste(response_name, "~ 1 + (1 | native_site/parent)"))
  # ---- 1 predictor ----
  f_native                <- as.formula(paste(response_name, "~ native_site + (1 | native_site/parent)"))
  f_transplant            <- as.formula(paste(response_name, "~ transplant_site + (1 | native_site/parent)"))
  f_time                  <- as.formula(paste(response_name, "~ timepoint + (1 | native_site/parent)"))
  f_sym                   <- as.formula(paste(response_name, "~ sym_dom + (1 | native_site/parent)"))
  # ---- 2 predictors, additive (all 6 pairs) ----
  f_site_add                  <- as.formula(paste(response_name, "~ native_site + transplant_site + (1 | native_site/parent)"))
  f_native_time               <- as.formula(paste(response_name, "~ native_site + timepoint + (1 | native_site/parent)"))                 # fixed: was missing "+"
  f_native_sym                <- as.formula(paste(response_name, "~ native_site + sym_dom + (1 | native_site/parent)"))
  f_transplant_time           <- as.formula(paste(response_name, "~ transplant_site + timepoint + (1 | native_site/parent)"))
  f_transplant_sym            <- as.formula(paste(response_name, "~ transplant_site + sym_dom + (1 | native_site/parent)"))
  f_sym_time                  <- as.formula(paste(response_name, "~ sym_dom + timepoint + (1 | native_site/parent)"))
  # ---- 2 predictors, interactive (all 6 pairs) ----
  f_site_int                  <- as.formula(paste(response_name, "~ native_site * transplant_site + (1 | native_site/parent)"))
  f_native_time_int           <- as.formula(paste(response_name, "~ native_site * timepoint + (1 | native_site/parent)"))
  f_native_sym_int            <- as.formula(paste(response_name, "~ native_site * sym_dom + (1 | native_site/parent)"))
  f_transplant_time_int       <- as.formula(paste(response_name, "~ transplant_site * timepoint + (1 | native_site/parent)"))
  f_transplant_sym_int        <- as.formula(paste(response_name, "~ transplant_site * sym_dom + (1 | native_site/parent)"))
  f_sym_time_int              <- as.formula(paste(response_name, "~ sym_dom * timepoint + (1 | native_site/parent)"))            # NEW (symmetry)
  # ---- 3 predictors, all additive (all 4 triples) ----
  f_site_add_time_add         <- as.formula(paste(response_name, "~ native_site + transplant_site + timepoint + (1 | native_site/parent)"))
  f_site_add_sym_add          <- as.formula(paste(response_name, "~ native_site + transplant_site + sym_dom + (1 | native_site/parent)"))
  f_sym_time_native           <- as.formula(paste(response_name, "~ sym_dom + timepoint + native_site + (1 | native_site/parent)"))
  f_sym_time_transplant       <- as.formula(paste(response_name, "~ sym_dom + timepoint + transplant_site + (1 | native_site/parent)"))
  # ---- 3 predictors, one pair interacts + third additive (all 6 pairs) ----
  f_site_int_time_add                         <- as.formula(paste(response_name, "~ native_site * transplant_site + timepoint + (1 | native_site/parent)"))
  f_site_int_sym_add                          <- as.formula(paste(response_name, "~ native_site * transplant_site + sym_dom + (1 | native_site/parent)"))
  f_transplant_time_add_sym_native_int        <- as.formula(paste(response_name, "~ transplant_site + timepoint + sym_dom * native_site + (1 | native_site/parent)"))
  f_transplant_sym_add_time_native_int        <- as.formula(paste(response_name, "~ transplant_site + sym_dom + timepoint * native_site + (1 | native_site/parent)"))
  f_native_time_add_sym_transplant_int        <- as.formula(paste(response_name, "~ native_site + timepoint + sym_dom * transplant_site + (1 | native_site/parent)"))
  f_native_sym_add_time_transplant_int        <- as.formula(paste(response_name, "~ native_site + sym_dom + timepoint * transplant_site + (1 | native_site/parent)"))
  # ---- 3 predictors, "hub" (one variable interacts with each of the
  # other two separately, those two don't interact with each other) ----
  f_site_add_sym_int                                 <- as.formula(paste(response_name, "~ (native_site + transplant_site) * sym_dom + (1 | native_site/parent)"))
  f_site_add_time_int                                <- as.formula(paste(response_name, "~ (native_site + transplant_site) * timepoint + (1 | native_site/parent)"))
  f_native_sym_add_time_int                          <- as.formula(paste(response_name, "~ (native_site + sym_dom) * timepoint + (1 | native_site/parent)"))               # NEW (hub)
  f_transplant_sym_add_time_int                      <- as.formula(paste(response_name, "~ (transplant_site + sym_dom) * timepoint + (1 | native_site/parent)"))          # NEW (hub)
  f_native_time_add_sym_int                          <- as.formula(paste(response_name, "~ (native_site + timepoint) * sym_dom + (1 | native_site/parent)"))               # NEW (hub)
  f_transplant_time_add_sym_int                      <- as.formula(paste(response_name, "~ (transplant_site + timepoint) * sym_dom + (1 | native_site/parent)")) 
  # ---- 3 predictors, fully interactive, 4th predictor absent (all 4 triples) ----
  f_site_int_sym_int                                    <- as.formula(paste(response_name, "~ native_site * transplant_site * sym_dom + (1 | native_site/parent)"))
  f_site_int_time_int                                   <- as.formula(paste(response_name, "~ native_site * transplant_site * timepoint + (1 | native_site/parent)"))
  f_native_sym_time_int                                 <- as.formula(paste(response_name, "~ native_site * sym_dom * timepoint + (1 | native_site/parent)"))             # NEW (symmetry)
  f_transplant_sym_time_int                             <- as.formula(paste(response_name, "~ transplant_site * sym_dom * timepoint + (1 | native_site/parent)"))        # NEW (symmetry)
  # ---- 4 predictors, fully additive ----
  f_sym_time_site_add                                   <- as.formula(paste(response_name, "~ sym_dom + timepoint + native_site + transplant_site + (1 | native_site/parent)"))
  # ---- 4 predictors, one pair interacts + other two additive ----
  f_site_int_sym_time_add                                     <- as.formula(paste(response_name, "~ native_site * transplant_site + sym_dom + timepoint + (1 | native_site/parent)"))     # fixed: stray ")"
  f_site_add_sym_time_int                                     <- as.formula(paste(response_name, "~ native_site + transplant_site + timepoint * sym_dom + (1 | native_site/parent)"))
  f_native_sym_int_transplant_time_add                        <- as.formula(paste(response_name, "~ native_site * sym_dom + timepoint + transplant_site + (1 | native_site/parent)"))     # fixed: stray ")"
  f_transplant_sym_int_native_time_add                        <- as.formula(paste(response_name, "~  transplant_site * sym_dom + timepoint + native_site + (1 | native_site/parent)"))
  f_native_time_int_transplant_sym_add                        <- as.formula(paste(response_name, "~ native_site * timepoint + transplant_site + sym_dom + (1 | native_site/parent)"))     # fixed: stray ")"
  f_transplant_time_int_native_sym_add                        <- as.formula(paste(response_name, "~ transplant_site * timepoint + sym_dom + native_site + (1 | native_site/parent)"))
  # ---- 4 predictors, one triple fully interacts + 4th additive (all 4 triples) ----
  f_site_sym_int_time_add                                    <- as.formula(paste(response_name, "~ native_site * transplant_site * sym_dom + timepoint + (1 | native_site/parent)"))
  f_site_time_int_sym_add                                    <- as.formula(paste(response_name, "~ native_site * transplant_site * timepoint + sym_dom + (1 | native_site/parent)"))
  f_native_sym_time_int_transplant_add                       <- as.formula(paste(response_name, "~ native_site * sym_dom * timepoint + transplant_site + (1 | native_site/parent)"))
  f_transplant_sym_time_int_native_add                       <- as.formula(paste(response_name, "~ transplant_site * sym_dom * timepoint + native_site + (1 | native_site/parent)"))
  # ---- 4 predictors, fully interactive ----
  f_site_sym_time_int                                       <- as.formula(paste(response_name, "~ native_site * transplant_site * timepoint * sym_dom + (1 | native_site/parent)"))  # fixed: was missing "+"

  models <- list(
    null                                  = fit_one(f_null, data),

    native                                = fit_one(f_native, data),
    transplant                            = fit_one(f_transplant, data),
    time                                  = fit_one(f_time, data),
    sym                                   = fit_one(f_sym, data),

    site_add                              = fit_one(f_site_add, data),
    native_time                           = fit_one(f_native_time, data),
    native_sym                            = fit_one(f_native_sym, data),
    transplant_time                       = fit_one(f_transplant_time, data),
    transplant_sym                        = fit_one(f_transplant_sym, data),
    sym_time                              = fit_one(f_sym_time, data),

    site_int                              = fit_one(f_site_int, data),
    native_time_int                       = fit_one(f_native_time_int, data),
    native_sym_int                        = fit_one(f_native_sym_int, data),
    transplant_time_int                   = fit_one(f_transplant_time_int, data),
    transplant_sym_int                    = fit_one(f_transplant_sym_int, data),
    sym_time_int                          = fit_one(f_sym_time_int, data),

    site_add_time_add                     = fit_one(f_site_add_time_add, data),
    site_add_sym_add                      = fit_one(f_site_add_sym_add, data),
    sym_time_native                       = fit_one(f_sym_time_native, data),
    sym_time_transplant                   = fit_one(f_sym_time_transplant, data),

    site_int_time_add                     = fit_one(f_site_int_time_add, data),
    site_int_sym_add                      = fit_one(f_site_int_sym_add, data),
    transplant_time_add_sym_native_int    = fit_one(f_transplant_time_add_sym_native_int, data),
    transplant_sym_add_time_native_int    = fit_one(f_transplant_sym_add_time_native_int, data),
    native_time_add_sym_transplant_int    = fit_one(f_native_time_add_sym_transplant_int, data),
    native_sym_add_time_transplant_int    = fit_one(f_native_sym_add_time_transplant_int, data),

    site_add_sym_int                      = fit_one(f_site_add_sym_int, data),
    site_add_time_int                     = fit_one(f_site_add_time_int, data),
    native_sym_add_time_int               = fit_one(f_native_sym_add_time_int, data),
    transplant_sym_add_time_int           = fit_one(f_transplant_sym_add_time_int, data),
    native_time_add_sym_int               = fit_one(f_native_time_add_sym_int, data),
    transplant_time_add_sym_int           = fit_one(f_transplant_time_add_sym_int, data),

    site_int_sym_int                      = fit_one(f_site_int_sym_int, data),
    site_int_time_int                     = fit_one(f_site_int_time_int, data),
    native_sym_time_int                   = fit_one(f_native_sym_time_int, data),
    transplant_sym_time_int               = fit_one(f_transplant_sym_time_int, data),

    sym_time_site_add                     = fit_one(f_sym_time_site_add, data),

    site_int_sym_time_add                 = fit_one(f_site_int_sym_time_add, data),
    site_add_sym_time_int                 = fit_one(f_site_add_sym_time_int, data),
    native_sym_int_transplant_time_add    = fit_one(f_native_sym_int_transplant_time_add, data),
    transplant_sym_int_native_time_add    = fit_one(f_transplant_sym_int_native_time_add, data),
    native_time_int_transplant_sym_add    = fit_one(f_native_time_int_transplant_sym_add, data),
    transplant_time_int_native_sym_add    = fit_one(f_transplant_time_int_native_sym_add, data),

    site_sym_int_time_add                 = fit_one(f_site_sym_int_time_add, data),
    site_time_int_sym_add                 = fit_one(f_site_time_int_sym_add, data),
    native_sym_time_int_transplant_add    = fit_one(f_native_sym_time_int_transplant_add, data),
    transplant_sym_time_int_native_add    = fit_one(f_transplant_sym_time_int_native_add, data),

    site_sym_time_int                     = fit_one(f_site_sym_time_int, data)
  )

  return(models)
}


####### Helper: build & fit the SAME 14 models, WITHOUT timepoint #########

fit_models_no_time <- function(response_name, data) {

  f_null               <- as.formula(paste(response_name, "~ 1 + (1 | native_site/parent)"))
  f_native             <- as.formula(paste(response_name, "~ native_site + (1 | native_site/parent)"))
  f_transplant         <- as.formula(paste(response_name, "~ transplant_site + (1 | native_site/parent)"))
  f_site_add           <- as.formula(paste(response_name, "~ native_site + transplant_site + (1 | native_site/parent)"))
  f_site_int           <- as.formula(paste(response_name, "~ native_site * transplant_site + (1 | native_site/parent)"))
  f_sym                <- as.formula(paste(response_name, "~ sym_dom + (1 | native_site/parent)"))
  f_native_sym         <- as.formula(paste(response_name, "~ native_site + sym_dom + (1 | native_site/parent)"))
  f_transplant_sym     <- as.formula(paste(response_name, "~ transplant_site + sym_dom + (1 | native_site/parent)"))
  f_site_add_sym_add   <- as.formula(paste(response_name, "~ native_site + transplant_site + sym_dom + (1 | native_site/parent)"))
  f_site_int_sym_add   <- as.formula(paste(response_name, "~ native_site * transplant_site + sym_dom + (1 | native_site/parent)"))
  f_native_sym_int     <- as.formula(paste(response_name, "~ native_site * sym_dom + (1 | native_site/parent)"))
  f_transplant_sym_int <- as.formula(paste(response_name, "~ transplant_site * sym_dom + (1 | native_site/parent)"))
  f_site_add_sym_int   <- as.formula(paste(response_name, "~ (native_site + transplant_site) * sym_dom + (1 | native_site/parent)"))
  f_site_int_sym_int   <- as.formula(paste(response_name, "~ native_site * transplant_site * sym_dom + (1 | native_site/parent)"))

  models <- list(
    null                 = fit_one(f_null, data),
    native               = fit_one(f_native, data),
    transplant           = fit_one(f_transplant, data),
    site_add             = fit_one(f_site_add, data),
    site_int             = fit_one(f_site_int, data),
    sym                  = fit_one(f_sym, data),
    native_sym           = fit_one(f_native_sym, data),
    transplant_sym       = fit_one(f_transplant_sym, data),
    site_add_sym_add     = fit_one(f_site_add_sym_add, data),
    site_int_sym_add     = fit_one(f_site_int_sym_add, data),
    native_sym_int       = fit_one(f_native_sym_int, data),
    transplant_sym_int   = fit_one(f_transplant_sym_int, data),
    site_add_sym_int     = fit_one(f_site_add_sym_int, data),
    site_int_sym_int     = fit_one(f_site_int_sym_int, data)
  )

  return(models)
}


####### Helper: compare a list of fitted models by AIC #########

compare_models_by_aic <- function(models_list, response_name, species_name) {

  good_models <- models_list[!sapply(models_list, is.null)]

  # If EVERY candidate model failed, return an empty (but properly
  # shaped) table instead of erroring.
  if (length(good_models) == 0) {
    return(data.frame(
      species = character(0), response = character(0), model = character(0),
      df = numeric(0), AIC = numeric(0), singular = logical(0),
      deltaAIC = numeric(0), weight = numeric(0)
    ))
  }

  model_names     <- names(good_models)
  aic_values      <- sapply(good_models, AIC)
  df_values       <- sapply(good_models, function(m) attr(logLik(m), "df"))
  singular_values <- sapply(good_models, isSingular)

  results <- data.frame(
    species  = species_name,
    response = response_name,
    model    = model_names,
    df       = df_values,
    AIC      = aic_values,
    singular = singular_values
  )

  results <- results[order(results$AIC), ]
  results$deltaAIC <- results$AIC - min(results$AIC)
  results$weight   <- exp(-0.5 * results$deltaAIC) / sum(exp(-0.5 * results$deltaAIC))

  rownames(results) <- NULL
  return(results)
}


####### Helper: simple jitter plots #########

plot_trait_with_time <- function(data, response_name, species_name) {
  ggplot(data, aes(x = timepoint, y = .data[[response_name]], color = sym_dom)) +
    geom_jitter(width = 0.15, alpha = 0.85, size = 2.5) +
    scale_color_manual(values = sym_dom_colors, name = "Dominant\nITS2 type") +
    facet_wrap(~ transplant_group, nrow = 1) +
    theme_bw(base_size = 15) +
    theme(axis.title.x = element_blank()) +
    labs(y = response_name, title = paste(response_name, "-", species_name))
}

plot_trait_no_time <- function(data, response_name, species_name) {
  ggplot(data, aes(x = transplant_group, y = .data[[response_name]], color = sym_dom)) +
    geom_jitter(width = 0.15, alpha = 0.85, size = 2.5) +
    scale_color_manual(values = sym_dom_colors, name = "Dominant\nITS2 type") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
    labs(y = response_name, title = paste(response_name, "- T0 to T12 change -", species_name))
}


########################################################################
####### PART 1: PLOTS ONLY (look at these before trusting any model) #######
########################################################################


for (response_name in response_vars_with_time) {
  for (species_name in names(species_data)) {

    this_data <- species_data[[species_name]] %>%
      filter(!is.na(.data[[response_name]]))

    if (nrow(this_data) == 0) {
      cat("No data at all for", response_name, "in", species_name, "- no plot made.\n")
      next
    }

    file_tag <- paste0(response_name, "_", gsub(" ", "_", species_name))
    p <- plot_trait_with_time(this_data, response_name, species_name)
    ggsave(here("Output", "Simple_Models", paste0("plot_", file_tag, ".pdf")),
           p, width = 12, height = 5)
  }
}

for (response_name in response_vars_no_time) {
  for (species_name in names(species_data)) {

    this_data <- species_data[[species_name]] %>%
      filter(timepoint == "T12") %>%
      filter(!is.na(.data[[response_name]]))

    if (nrow(this_data) == 0) {
      cat("No data at all for", response_name, "in", species_name, "- no plot made.\n")
      next
    }

    file_tag <- paste0(response_name, "_", gsub(" ", "_", species_name))
    p <- plot_trait_no_time(this_data, response_name, species_name)
    ggsave(here("Output", "Simple_Models", paste0("plot_", file_tag, ".pdf")),
           p, width = 8, height = 5)
  }
}



########################################################################
####### PART 2: FIT + CHECK + AIC COMPARE #######
########################################################################
# Everything gets stored in these two master lists so you can pull any
# individual model or table back out afterward without re-running
# anything, e.g.:
#   all_fitted_models[["PNET"]][["Porites sp."]]$site_int
#   all_aic_tables[["PNET"]][["Porites sp."]]

all_fitted_models <- list()
all_aic_tables    <- list()

# ---- traits WITH timepoint ----
for (response_name in response_vars_with_time) {

  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models[[response_name]] <- list()
  all_aic_tables[[response_name]]    <- list()

  for (species_name in names(species_data)) {

    cat("\n---", species_name, "---\n")

    this_data <- species_data[[species_name]] %>%
      filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
      droplevels()

    cat("Rows of data after filtering:", nrow(this_data), "\n")

    if (nrow(this_data) < 10) {
      cat("Not enough data - skipping.\n")
      next
    }

    # ---- FIT all 43 candidate models ----
    models_list <- fit_models_with_time(response_name, this_data)
    all_fitted_models[[response_name]][[species_name]] <- models_list

    n_fit    <- sum(!sapply(models_list, is.null))
    n_failed <- sum(sapply(models_list, is.null))
    cat(n_fit, "of 43 models fit successfully;", n_failed, "failed.\n")

    # ---- AIC COMPARE ----
    aic_table <- compare_models_by_aic(models_list, response_name, species_name)
    all_aic_tables[[response_name]][[species_name]] <- aic_table

    if (nrow(aic_table) == 0) {
      cat("None of the 43 models fit for", response_name, "in", species_name, "- see [ERROR] messages above.\n")
      next
    }

    print(aic_table)

    file_tag <- paste0(response_name, "_", gsub(" ", "_", species_name))
    write_csv(aic_table, here("Output", "Simple_Models", paste0("aic_table_", file_tag, ".csv")))

    # ---- CHECK MODEL on the best-by-AIC model ----
    best_model_name <- aic_table$model[1]
    best_model      <- models_list[[best_model_name]]
    cat("Best model:", best_model_name, "\n")

    best_model_reml <- tryCatch(
      lmer(formula(best_model), data = this_data, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) {
        cat("Could not refit best model with REML:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(best_model_reml)) {
      tryCatch({
        pdf(here("Output", "Simple_Models", paste0("checkmodel_", file_tag, ".pdf")), width = 12, height = 10)
        print(check_model(best_model_reml))
      }, error = function(e) {
        cat("check_model() failed:", conditionMessage(e), "\n")
      }, finally = {
        dev.off()
      })
    }
  }
}

# ---- traits WITHOUT timepoint (T12 only) ----
for (response_name in response_vars_no_time) {

  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models[[response_name]] <- list()
  all_aic_tables[[response_name]]    <- list()

  for (species_name in names(species_data)) {

    cat("\n---", species_name, "---\n")

    this_data <- species_data[[species_name]] %>%
      filter(timepoint == "T12") %>%
      filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
      droplevels()

    cat("Rows of data after filtering:", nrow(this_data), "\n")

    if (nrow(this_data) < 10) {
      cat("Not enough data - skipping.\n")
      next
    }

    # ---- FIT all 14 candidate models ----
    models_list <- fit_models_no_time(response_name, this_data)
    all_fitted_models[[response_name]][[species_name]] <- models_list

    n_fit    <- sum(!sapply(models_list, is.null))
    n_failed <- sum(sapply(models_list, is.null))
    cat(n_fit, "of 14 models fit successfully;", n_failed, "failed.\n")

    # ---- AIC COMPARE ----
    aic_table <- compare_models_by_aic(models_list, response_name, species_name)
    all_aic_tables[[response_name]][[species_name]] <- aic_table

    if (nrow(aic_table) == 0) {
      cat("None of the 14 models fit for", response_name, "in", species_name, "- see [ERROR] messages above.\n")
      next
    }

    print(aic_table)

    file_tag <- paste0(response_name, "_", gsub(" ", "_", species_name))
    write_csv(aic_table, here("Output", "Simple_Models", paste0("aic_table_", file_tag, ".csv")))

    # ---- CHECK MODEL on the best-by-AIC model ----
    best_model_name <- aic_table$model[1]
    best_model      <- models_list[[best_model_name]]
    cat("Best model:", best_model_name, "\n")

    best_model_reml <- tryCatch(
      lmer(formula(best_model), data = this_data, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa")),
      error = function(e) {
        cat("Could not refit best model with REML:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(best_model_reml)) {
      tryCatch({
        pdf(here("Output", "Simple_Models", paste0("checkmodel_", file_tag, ".pdf")), width = 12, height = 10)
        print(check_model(best_model_reml))
      }, error = function(e) {
        cat("check_model() failed:", conditionMessage(e), "\n")
      }, finally = {
        dev.off()
      })
    }
  }
}


########################################################################
####### PART 3: WHICH PREDICTORS WIN - SUMMARY TABLE #######
########################################################################
# No refitting here -- just pulling the top (best-AIC) row out of each
# table we already saved in Part 2.

all_response_vars <- c(response_vars_with_time, response_vars_no_time)

summary_rows <- list()

for (response_name in all_response_vars) {
  for (species_name in names(species_data)) {

    aic_table <- all_aic_tables[[response_name]][[species_name]]

    if (is.null(aic_table) || nrow(aic_table) == 0) next

    top_row <- aic_table[1, ]

    # Pull the actual fitted model back out of Part 2's storage and grab
    # its formula, rather than looking up a plain-English description.
    winning_model <- all_fitted_models[[response_name]][[species_name]][[top_row$model]]
    top_row$model_formula <- paste(deparse(formula(winning_model)), collapse = " ")

    summary_rows[[length(summary_rows) + 1]] <- top_row
  }
}

which_predictors_win <- do.call(rbind, summary_rows)

print(which_predictors_win)

write_csv(which_predictors_win, here("Output", "Simple_Models", "which_predictors_win_summary.csv"))
