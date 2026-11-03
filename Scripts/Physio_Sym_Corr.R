####### Physiology & symbiont ID correlations for Curaçao RTE 2026 #######
# Maya Powell
# Adapted with Claude: Aug 19 2026
#
# GOAL: ONE unified candidate model set per trait x species that treats
# native_site, transplant_site, sym_dom_b4, and timepoint all as first-class
# predictors (rather than a two-stage cascade). The catalog systematically
# covers:
#   - site structure: none / native_site only / transplant_site only /
#     native+transplant additive / native*transplant interactive
#   - symbiont: absent / present (and, when site is also present, whether
#     site and symbiont combine additively or interactively)
#   - time: additive (main effect only) / interactive (the whole predictor
#     structure above is allowed to change across timepoints)
# giving up to 27 candidate models (26 + null), all compared by AIC in a
# single contest per trait x species. Saturated combinations (symbiont
# crossed with site and/or time) are only attempted when there's enough
# replication to support them (see `guard`).
#
# Section 0: libraries, colors, data read-in, SPECIES SPLIT
# Section 1: exploratory jitter plots
# Section 2: shared helpers
# Section 3: the unified model catalog
# Section 4: formula builder + replication guard
# Section 5: fit & AIC-compare the full candidate set per trait x species
# Section 6: master "which predictor(s) win" table, with a readable label
#            auto-generated from the catalog metadata
# Section 7: refit winners for inference (Anova, R2, VIF)
# Section 8: compact letter displays (CLD) for the site term, time-aware
# Section 9: site plots


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
library(performance)
library(purrr)
library(rlang)
library(tibble)
library(psych)
library(corrplot)
library(ggpubr)
library(lme4)
library(multcompView)   # for compact letter displays (CLD) from pairwise p-values

####### Color palettes #########
sym_dom_b4_colors <- c("A4" = "#ffaabb", "C1" = "#eedd88", "C3" = "#77aadd",
                     "C_other" = "#99ddff", "D1" = "#994455")

site_colors <- c("Bay Initial" = "coral4", "Bay Native" = "coral4",
                  "Bay Transplant" = "salmon2", "Reef Transplant" = "lightblue3",
                  "Reef Initial" = "dodgerblue4", "Reef Native" = "dodgerblue4")

###### Read in data ########
sym_physio <- read_csv(here("Data", "sym_physio.csv"),
                        locale = locale(encoding = "latin1"))

###### Clean & prep data ########

sym_physio <- sym_physio %>%
  mutate(
    coral_species   = factor(coral_species),
    native_site     = factor(native_site, levels = c("Reef", "Bay")),
    transplant_site = factor(transplant_site, levels = c("Reef", "Bay")),
    sym_dom_b4         = factor(sym_dom_b4),
    timepoint       = factor(timepoint, levels = c("T0", "T4", "T12")),
    parent          = factor(parent),
    native_transplant = factor(native_transplant, levels = c("Initial", "Native", "Transplant")),
    transplant_group = factor(
      transplant_group,
      levels = c("Reef Initial", "Bay Initial",
                 "Reef Native", "Reef Transplant",
                 "Bay Native", "Bay Transplant")
    )
  )

####### SPECIES SPLIT ########
ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()

species_data <- list(
  "Siderastrea siderea" = ss_data,
  "Porites sp."          = pp_data
)

response_vars <- c("chla", "chla_sym", "dry_weight", "sym_density",
                    "calcification", "pct_change_dw", "pct_change_sa")

response_labels <- list(
  chla          = "Chlorophyll a (µg/cm2)",
  chla_sym      = "Chlorophyll a per symbiont (pg/cell)",
  dry_weight    = "Tissue biomass (mg/cm2)",
  sym_density   = "Symbiont density (x10^6 cells/cm2)",
  calcification = "Calcification (mg/cm2/d)",
  pct_change_dw = "% change tissue biomass (T0-T12)",
  pct_change_sa = "% change surface area (T0-T12)"
)

dir.create(here("Output", "Physiology"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("Output", "Tables"), recursive = TRUE, showWarnings = FALSE)

####### 1. Exploratory visualization (raw jitter, all data) ########

plot_physio_by_sym <- function(data, response) {
  ylab_txt <- response_labels[[response]] %||% response
  ggplot(data, aes(x = transplant_group, y = .data[[response]], color = sym_dom_b4)) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
    scale_color_manual(values = sym_dom_b4_colors) +
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
      here("Output", "Physiology", paste0(sp, "_", resp, "_by_symdom_jitter.pdf")),
      p, height = 5, width = 12
    )
  })
})

####### 2. Shared helpers ########

clean_trait_data <- function(data, response) {
  data %>%
    filter(
      !is.na(.data[[response]]), !is.na(sym_dom_b4),
      !is.na(native_site), !is.na(transplant_site)
    ) %>%
    droplevels()
}

fit_one <- function(f, d) {
  tryCatch(
    suppressWarnings(suppressMessages(
      lmer(f, data = d, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    )),
    error = function(e) NULL
  )
}

####### 3. The unified model catalog ########
# One row per candidate model. `core_term` is the fixed-effect formula
# fragment excluding timepoint and the random effect; `time_mode` says
# whether timepoint is added on as a main effect ("additive") or allowed
# to interact with everything in core_term ("interactive").
#
# has_native / has_transplant / has_sym flag which raw predictors are in
# core_term. site_internal_interactive flags whether native_site and
# transplant_site interact WITH EACH OTHER (vs. additive/absent).
# site_sym_int flags whether the site term and sym_dom_b4 interact with each
# other (vs. additive/absent).

model_catalog <- tribble(
  ~model,                    ~core_term,                                    ~time_mode,     ~has_native, ~has_transplant, ~has_sym, ~site_internal_interactive, ~site_sym_int,
  "null",                    "1",                                           "additive",     FALSE, FALSE, FALSE, FALSE, FALSE,

  "native",                  "native_site",                                 "additive",     TRUE,  FALSE, FALSE, FALSE, FALSE,
  "native_tp",               "native_site",                                 "interactive",  TRUE,  FALSE, FALSE, FALSE, FALSE,

  "transplant",              "transplant_site",                             "additive",     FALSE, TRUE,  FALSE, FALSE, FALSE,
  "transplant_tp",           "transplant_site",                             "interactive",  FALSE, TRUE,  FALSE, FALSE, FALSE,

  "site_add",                "native_site + transplant_site",               "additive",     TRUE,  TRUE,  FALSE, FALSE, FALSE,
  "site_add_tp",             "native_site + transplant_site",               "interactive",  TRUE,  TRUE,  FALSE, FALSE, FALSE,

  "site_int",                "native_site * transplant_site",               "additive",     TRUE,  TRUE,  FALSE, TRUE,  FALSE,
  "site_int_tp",             "native_site * transplant_site",               "interactive",  TRUE,  TRUE,  FALSE, TRUE,  FALSE,

  "sym",                     "sym_dom_b4",                                     "additive",     FALSE, FALSE, TRUE,  FALSE, FALSE,
  "sym_tp",                  "sym_dom_b4",                                     "interactive",  FALSE, FALSE, TRUE,  FALSE, FALSE,

  "native_sym",              "native_site + sym_dom_b4",                       "additive",     TRUE,  FALSE, TRUE,  FALSE, FALSE,
  "native_sym_tp",           "native_site + sym_dom_b4",                       "interactive",  TRUE,  FALSE, TRUE,  FALSE, FALSE,

  "transplant_sym",          "transplant_site + sym_dom_b4",                   "additive",     FALSE, TRUE,  TRUE,  FALSE, FALSE,
  "transplant_sym_tp",       "transplant_site + sym_dom_b4",                   "interactive",  FALSE, TRUE,  TRUE,  FALSE, FALSE,

  "site_add_sym_add",        "native_site + transplant_site + sym_dom_b4",     "additive",     TRUE,  TRUE,  TRUE,  FALSE, FALSE,
  "site_add_sym_add_tp",     "native_site + transplant_site + sym_dom_b4",     "interactive",  TRUE,  TRUE,  TRUE,  FALSE, FALSE,

  "site_int_sym_add",        "native_site * transplant_site + sym_dom_b4",     "additive",     TRUE,  TRUE,  TRUE,  TRUE,  FALSE,
  "site_int_sym_add_tp",     "native_site * transplant_site + sym_dom_b4",     "interactive",  TRUE,  TRUE,  TRUE,  TRUE,  FALSE,

  "native_sym_int",          "native_site * sym_dom_b4",                       "additive",     TRUE,  FALSE, TRUE,  FALSE, TRUE,
  "native_sym_int_tp",       "native_site * sym_dom_b4",                       "interactive",  TRUE,  FALSE, TRUE,  FALSE, TRUE,

  "transplant_sym_int",      "transplant_site * sym_dom_b4",                   "additive",     FALSE, TRUE,  TRUE,  FALSE, TRUE,
  "transplant_sym_int_tp",   "transplant_site * sym_dom_b4",                   "interactive",  FALSE, TRUE,  TRUE,  FALSE, TRUE,

  "site_add_sym_int",        "(native_site + transplant_site) * sym_dom_b4",   "additive",     TRUE,  TRUE,  TRUE,  FALSE, TRUE,
  "site_add_sym_int_tp",     "(native_site + transplant_site) * sym_dom_b4",   "interactive",  TRUE,  TRUE,  TRUE,  FALSE, TRUE,

  "site_int_sym_int",        "native_site * transplant_site * sym_dom_b4",     "additive",     TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
  "site_int_sym_int_tp",     "native_site * transplant_site * sym_dom_b4",     "interactive",  TRUE,  TRUE,  TRUE,  TRUE,  TRUE
) %>%
  mutate(
    site_structure = case_when(
      !has_native & !has_transplant ~ "none",
      has_native  & !has_transplant ~ "native",
      !has_native & has_transplant  ~ "transplant",
      has_native  & has_transplant & site_internal_interactive  ~ "interactive",
      has_native  & has_transplant & !site_internal_interactive ~ "additive"
    ),
    # Only fit the pricier combinations (symbiont crossed with site and/or
    # time) when there's enough replication to support them.
    guard = has_sym & (site_sym_int | time_mode == "interactive")
  )

####### 4. Formula builder + replication guard ########

build_formula <- function(response, core_term, time_mode, include_timepoint) {
  if (!include_timepoint) {
    return(as.formula(paste0(response, " ~ ", core_term, " + (1|parent)")))
  }
  if (time_mode == "additive") {
    as.formula(paste0(response, " ~ ", core_term, " + timepoint + (1|parent)"))
  } else {
    as.formula(paste0(response, " ~ (", core_term, ") * timepoint + (1|parent)"))
  }
}

guard_vars <- function(has_native, has_transplant, site_sym_int, time_mode) {
  vars <- "sym_dom_b4"
  if (site_sym_int) {
    if (has_native)     vars <- c(vars, "native_site")
    if (has_transplant) vars <- c(vars, "transplant_site")
  }
  if (time_mode == "interactive") vars <- c(vars, "timepoint")
  unique(vars)
}

enough_replication_generic <- function(d, vars, min_n = 2) {
  if (length(vars) == 0) return(TRUE)
  cell_n <- d %>% count(across(all_of(vars))) %>% pull(n)
  length(cell_n) > 0 && min(cell_n) >= min_n
}

####### 5. Fit & AIC-compare the full candidate set ########

fit_model_set <- function(data, response) {
  d <- clean_trait_data(data, response)
  if (nrow(d) < 10) return(list(data = d, fits = list(), catalog = model_catalog[0, ]))

  include_tp <- n_distinct(d$timepoint) > 1

  rows <- model_catalog %>% filter(time_mode == "additive" | include_tp)

  keep <- map_lgl(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    if (!r$guard) return(TRUE)
    enough_replication_generic(d, guard_vars(r$has_native, r$has_transplant, r$site_sym_int, r$time_mode))
  })
  rows <- rows[keep, ]

  fits <- pmap(rows, function(model, core_term, time_mode, ...) {
    f <- build_formula(response, core_term, time_mode, include_tp)
    fit_one(f, d)
  })
  names(fits) <- rows$model
  fits <- compact(fits)

  list(data = d, fits = fits, catalog = rows)
}

aic_table <- function(fit_list, response, species) {
  fits <- fit_list$fits
  if (length(fits) == 0) return(NULL)
  tibble(
    response = response, species = species, model = names(fits), n_obs = nrow(fit_list$data),
    df       = map_dbl(fits, ~ attr(logLik(.x), "df")),
    AIC      = map_dbl(fits, AIC),
    singular = map_lgl(fits, isSingular)
  ) %>%
    left_join(
      fit_list$catalog %>%
        select(model, core_term, time_mode, site_structure, has_sym, site_sym_int),
      by = "model"
    ) %>%
    mutate(
      deltaAIC = AIC - min(AIC),
      weight   = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
    ) %>%
    arrange(AIC)
}

all_fits <- imap(species_data, function(d, sp) {
  set_names(map(response_vars, ~ fit_model_set(d, .x)), response_vars)
})

model_comparison_all <- imap_dfr(species_data, function(d, sp) {
  map_dfr(response_vars, function(resp) {
    aic_table(all_fits[[sp]][[resp]], resp, sp)
  })
})

write_csv(model_comparison_all, here("Output", "Tables", "symbiont_physio_model_comparison_full_sym_b4.csv"))

####### 6. Master "which predictor(s) win" table ########

describe_model <- function(model, site_structure, has_sym, site_sym_int, time_mode) {
  if (model == "null") return("No predictor outperforms time alone")

  site_txt <- switch(site_structure,
    none        = NA_character_,
    native      = "Native site (origin)",
    transplant  = "Transplant site (destination)",
    additive    = "Native + Transplant site (additive, independent effects)",
    interactive = "Native x Transplant site (true reciprocal-transplant effect)"
  )
  sym_txt <- if (has_sym) "Symbiont identity" else NA_character_

  combo_txt <- if (!is.na(site_txt) && !is.na(sym_txt)) {
    paste0(site_txt, if (site_sym_int) " x " else " + ", sym_txt)
  } else {
    paste(na.omit(c(site_txt, sym_txt)), collapse = "")
  }

  time_txt <- if (time_mode == "interactive") " x Time (effect changes across timepoints)" else " (constant across time)"
  paste0(combo_txt, time_txt)
}

best_models <- model_comparison_all %>%
  group_by(response, species) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    predictor_support = describe_model(model, site_structure, has_sym, site_sym_int, time_mode),
    time_varying = time_mode == "interactive" & model != "null",
    meaningfully_better_than_null = deltaAIC == 0 & model != "null"
  ) %>%
  ungroup()

print(best_models %>%
        select(species, response, model, predictor_support, site_structure, has_sym, time_varying, AIC, deltaAIC, weight))

write_csv(best_models, here("Output", "Tables", "symbiont_physio_best_models_summary_sym_b4.csv"))

####### 7. Refit winning models with REML for inference ########

refit_best <- function(data, response, core_term, time_mode) {
  d <- clean_trait_data(data, response)
  include_tp <- n_distinct(d$timepoint) > 1
  f <- build_formula(response, core_term, time_mode, include_tp)
  tryCatch(
    lmer(f, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
}

inference_table <- best_models %>% filter(model != "null")

inference_results <- pmap(inference_table, function(response, species, model, core_term, time_mode,
                                                       site_structure, has_sym, site_sym_int, ...) {
  d   <- species_data[[species]]
  mod <- refit_best(d, response, core_term, time_mode)
  if (is.null(mod)) return(NULL)

  aov_tab <- tryCatch(car::Anova(mod), error = function(e) NULL)   # Type II Wald chi-sq
  r2_tab  <- tryCatch(performance::r2(mod), error = function(e) NULL)
  has_site <- site_structure != "none"
  vif_tab <- if (has_site && has_sym) {
    tryCatch(car::vif(mod), error = function(e) NULL)
  } else NULL

  list(response = response, species = species, model_name = model,
       site_structure = site_structure, time_mode = time_mode,
       model_obj = mod, anova = aov_tab, r2 = r2_tab, vif = vif_tab)
})
names(inference_results) <- paste(inference_table$species, inference_table$response, sep = " | ")
inference_results <- compact(inference_results)

walk(inference_results, function(res) {
  cat("\n==========================================================\n")
  cat(res$species, "-", res$response, "| Best model:", res$model_name, "\n")
  cat("==========================================================\n")
  print(res$anova)
  cat("\n")
  print(res$r2)
  if (!is.null(res$vif)) {
    cat("\nVariance inflation (collinearity check, site vs. symbiont terms):\n")
    print(res$vif)
  }
})

####### 8. Compact letter displays (CLD) for the site term ########
# For plotting, always display groups using the human-readable
# `transplant_group` factor (Reef Native / Reef Transplant / Bay Native /
# Bay Transplant / ...), refitting the winning model's *time structure*
# with transplant_group substituted for whatever site term won (native
# only / transplant only / additive / interactive) so labels/colors line
# up with `site_colors`. The subtitle still reports the TRUE winning site
# structure from Section 6/7 so the plot never overstates what the model
# actually found.

get_cld <- function(model, group_var, by_var = NULL) {
  tryCatch({
    if (is.null(by_var)) {
      emm <- emmeans(model, as.formula(paste("~", group_var)))
      pw  <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
        mutate(contrast_clean = str_replace_all(contrast, " - ", "-"))
      pvec <- setNames(pw$p.value, pw$contrast_clean)
      letters_out <- multcompLetters(pvec)$Letters
      tibble(group = names(letters_out), letter = unname(letters_out))
    } else {
      emm <- emmeans(model, as.formula(paste("~", group_var, "|", by_var)))
      pw  <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
        mutate(contrast_clean = str_replace_all(contrast, " - ", "-"))
      pw %>%
        group_by(.data[[by_var]]) %>%
        group_modify(~ {
          pvec <- setNames(.x$p.value, .x$contrast_clean)
          letters_out <- multcompLetters(pvec)$Letters
          tibble(group = names(letters_out), letter = unname(letters_out))
        }) %>%
        ungroup()
    }
  }, error = function(e) NULL)
}

term_is_significant <- function(anova_tab, pattern, alpha = 0.05) {
  if (is.null(anova_tab)) return(FALSE)
  rows <- rownames(anova_tab)
  hit  <- grepl(pattern, rows, fixed = TRUE)
  any(anova_tab[["Pr(>Chisq)"]][hit] < alpha, na.rm = TRUE)
}

plot_ready <- pmap(inference_table, function(response, species, model, core_term, time_mode,
                                              site_structure, ...) {
  d <- species_data[[species]]
  mod_plot <- refit_best(d, response, "transplant_group", time_mode)
  if (is.null(mod_plot)) return(NULL)

  aov_plot <- tryCatch(car::Anova(mod_plot), error = function(e) NULL)

  site_present <- site_structure != "none"
  site_sig <- site_present && term_is_significant(aov_plot, "transplant_group")

  cld_site <- if (site_sig) {
    if (time_mode == "interactive") get_cld(mod_plot, "transplant_group", "timepoint")
    else get_cld(mod_plot, "transplant_group")
  } else NULL

  list(
    response = response, species = species, model_name = model,
    site_structure = site_structure, time_mode = time_mode,
    model_obj = mod_plot, anova = aov_plot,
    site_sig = site_sig, site_time_varying = time_mode == "interactive",
    cld_site = cld_site
  )
})
names(plot_ready) <- paste(inference_table$species, inference_table$response, sep = " | ")
plot_ready <- compact(plot_ready)

####### 9. Site plots ########

dir.create(here("Output", "Physiology", "Site"), recursive = TRUE, showWarnings = FALSE)

subtitle_for <- function(res) {
  paste0("Best model: ", res$model_name, " (site structure: ", res$site_structure, ")",
         " | ", res$species,
         if (!is.null(res$anova)) {
           sig <- rownames(res$anova)[res$anova[["Pr(>Chisq)"]] < 0.05]
           if (length(sig) > 0) paste0(" | sig: ", paste(sig, collapse = ", ")) else " | no sig. terms"
         } else "")
}

drop_levels_per_facet <- function(data, facet_var, level_var) {
  data %>%
    group_by(.data[[facet_var]]) %>%
    group_modify(~ mutate(.x, "{level_var}" := fct_drop(.data[[level_var]]))) %>%
    ungroup()
}

join_cld <- function(lab_pos, cld, group_col) {
  if ("timepoint" %in% names(cld)) {
    lab_pos %>% left_join(cld, by = c("timepoint", setNames("group", group_col)))
  } else {
    lab_pos %>% left_join(cld, by = setNames("group", group_col))
  }
}

plot_site_trait <- function(res) {
  d <- clean_trait_data(species_data[[res$species]], res$response) %>%
    drop_levels_per_facet("timepoint", "transplant_group")

  ylab_txt <- response_labels[[res$response]] %||% res$response

  p <- ggplot(d, aes(x = transplant_group, y = .data[[res$response]])) +
    geom_boxplot(aes(fill = transplant_group), alpha = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = sym_dom_b4), width = 0.15, alpha = 0.85, size = 2) +
    scale_fill_manual(values = site_colors, guide = "none") +
    scale_color_manual(values = sym_dom_b4_colors, name = "Dominant\nITS2 type") +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.subtitle = element_text(size = 10, face = "italic")
    ) +
    labs(y = ylab_txt, title = res$response, subtitle = subtitle_for(res)) +
    facet_wrap(~ timepoint, scales = "free")

  if (!is.null(res$cld_site)) {
    lab_pos <- d %>%
      group_by(timepoint, transplant_group) %>%
      summarise(y = max(.data[[res$response]], na.rm = TRUE) * 1.08, .groups = "drop") %>%
      join_cld(res$cld_site, "transplant_group")

    p <- p + geom_text(data = lab_pos, aes(x = transplant_group, y = y, label = letter),
                        inherit.aes = FALSE, fontface = "bold", size = 5)
  }
  p
}

walk(plot_ready, function(res) {
  p_site <- plot_site_trait(res)
  ggsave(
    here("Output", "Physiology", "Site", paste0(res$species, "_", res$response, "_site.pdf")),
    p_site, height = 5.5, width = 10
  )
})

# View one interactively, e.g.:
# plot_site_trait(plot_ready[["Siderastrea siderea | dry_weight"]])
