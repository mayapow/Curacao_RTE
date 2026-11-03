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
emm_example("cld-multcomp")
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
    transplant_group = factor(
      transplant_group,
      levels = c("Reef Initial", "Bay Initial",
                 "Reef Native", "Reef Transplant",
                 "Bay Native", "Bay Transplant")))

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


#### Modeling setup functions for data both across and within timepoints ####

#clean up trait data to make sure all variables are good and have values
clean_trait_data <- function(data, response) {
  data %>%
    filter(
      !is.na(.data[[response]]), !is.na(sym_dom),
      !is.na(native_site), !is.na(native_transplant)
    ) %>%
    droplevels()
}

## fit helper with optimizer
fit_one <- function(f, d) {
  tryCatch(
    suppressWarnings(suppressMessages(
      lmer(f, data = d, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    )),
    error = function(e) NULL
  )
}

##write out all models for all variables with time
fit_models_with_time <- function(response_name, data) {
  # ---- 0 predictors ----
  f_null                 <- as.formula(paste(response_name, "~ 1 + (1 | native_site:parent)"))
  # ---- 1 predictor ----
  f_native                <- as.formula(paste(response_name, "~ native_site + (1 | native_site:parent)"))
  f_transplant            <- as.formula(paste(response_name, "~ native_transplant + (1 | native_site:parent)"))
  f_time                  <- as.formula(paste(response_name, "~ timepoint + (1 | native_site:parent)"))
  f_sym                   <- as.formula(paste(response_name, "~ sym_dom + (1 | native_site:parent)"))
  # ---- 2 predictors, additive (all 6 pairs) ----
  f_site_add                  <- as.formula(paste(response_name, "~ native_site + native_transplant + (1 | native_site:parent)"))
  f_native_time               <- as.formula(paste(response_name, "~ native_site + timepoint + (1 | native_site:parent)"))                 # fixed: was missing "+"
  f_native_sym                <- as.formula(paste(response_name, "~ native_site + sym_dom + (1 | native_site:parent)"))
  f_transplant_time           <- as.formula(paste(response_name, "~ native_transplant + timepoint + (1 | native_site:parent)"))
  f_transplant_sym            <- as.formula(paste(response_name, "~ native_transplant + sym_dom + (1 | native_site:parent)"))
  f_sym_time                  <- as.formula(paste(response_name, "~ sym_dom + timepoint + (1 | native_site:parent)"))
  # ---- 2 predictors, interactive (all 6 pairs) ----
  f_site_int                  <- as.formula(paste(response_name, "~ native_site * native_transplant + (1 | native_site:parent)"))
  f_native_time_int           <- as.formula(paste(response_name, "~ native_site * timepoint + (1 | native_site:parent)"))
  f_native_sym_int            <- as.formula(paste(response_name, "~ native_site * sym_dom + (1 | native_site:parent)"))
  f_transplant_time_int       <- as.formula(paste(response_name, "~ native_transplant * timepoint + (1 | native_site:parent)"))
  f_transplant_sym_int        <- as.formula(paste(response_name, "~ native_transplant * sym_dom + (1 | native_site:parent)"))
  f_sym_time_int              <- as.formula(paste(response_name, "~ sym_dom * timepoint + (1 | native_site:parent)"))            # NEW (symmetry)
  # ---- 3 predictors, all additive (all 4 triples) ----
  f_site_add_time_add         <- as.formula(paste(response_name, "~ native_site + native_transplant + timepoint + (1 | native_site:parent)"))
  f_site_add_sym_add          <- as.formula(paste(response_name, "~ native_site + native_transplant + sym_dom + (1 | native_site:parent)"))
  f_sym_time_native           <- as.formula(paste(response_name, "~ sym_dom + timepoint + native_site + (1 | native_site:parent)"))
  f_sym_time_transplant       <- as.formula(paste(response_name, "~ sym_dom + timepoint + native_transplant + (1 | native_site:parent)"))
  # ---- 3 predictors, one pair interacts + third additive ----
  f_site_int_time_add                         <- as.formula(paste(response_name, "~ native_site * native_transplant + timepoint + (1 | native_site:parent)"))
  f_site_int_sym_add                          <- as.formula(paste(response_name, "~ native_site * native_transplant + sym_dom + (1 | native_site:parent)"))
  f_native_sym_int_time_add                   <- as.formula(paste(response_name, "~ native_site * sym_dom + timepoint + (1 | native_site:parent)"))
  f_native_sym_int_transplant_add             <- as.formula(paste(response_name, "~ native_site * sym_dom + native_transplant + (1 | native_site:parent)"))
  f_native_time_int_sym_add                   <- as.formula(paste(response_name, "~ native_site * timepoint + sym_dom + (1 | native_site:parent)"))
  f_native_time_int_transplant_add            <- as.formula(paste(response_name, "~ native_site * timepoint + transplant + (1 | native_site:parent)"))
  f_transplant_sym_int_time_add               <- as.formula(paste(response_name, "~ native_transplant * sym_dom + timepoint + (1 | native_site:parent)"))
  f_transplant_sym_int_native_add             <- as.formula(paste(response_name, "~ native_transplant * sym_dom + native_site + (1 | native_site:parent)"))
  f_transplant_time_int_sym_add               <- as.formula(paste(response_name, "~ native_transplant * timepoint + sym_dom +(1 | native_site:parent)"))
  f_transplant_time_int_native_add            <- as.formula(paste(response_name, "~ native_transplant * timepoint + native_site + (1 | native_site:parent)"))
  f_sym_time_int_native_add                   <- as.formula(paste(response_name, "~ sym_dom * timepoint + native_site + (1 | native_site:parent)"))
  f_sym_time_int_transplant_add               <- as.formula(paste(response_name, "~ sym_dom * timepoint + native_transplant + (1 | native_site:parent)"))
  # ---- 3 predictors, "hub" (one variable interacts with each of the
  # other two separately, those two don't interact with each other) ----
  f_site_add_sym_int                                 <- as.formula(paste(response_name, "~ (native_site + native_transplant) * sym_dom + (1 | native_site:parent)"))
  f_site_add_time_int                                <- as.formula(paste(response_name, "~ (native_site + native_transplant) * timepoint + (1 | native_site:parent)"))
  f_native_sym_add_time_int                          <- as.formula(paste(response_name, "~ (native_site + sym_dom) * timepoint + (1 | native_site:parent)"))
  f_native_sym_add_transplant_int                    <- as.formula(paste(response_name, "~ (native_site + sym_dom) * native_transplant + (1 | native_site:parent)"))
  f_native_time_add_sym_int                          <- as.formula(paste(response_name, "~ (native_site + timepoint) * sym_dom + (1 | native_site:parent)"))  
  f_native_time_add_transplant_int                   <- as.formula(paste(response_name, "~ (native_site + timepoint) * native_transplant + (1 | native_site:parent)"))  
  f_transplant_sym_add_time_int                      <- as.formula(paste(response_name, "~ (native_transplant + sym_dom) * timepoint + (1 | native_site:parent)"))          
  f_transplant_sym_add_native_int                    <- as.formula(paste(response_name, "~ (native_transplant + sym_dom) * native_site + (1 | native_site:parent)"))          
  f_transplant_time_add_sym_int                      <- as.formula(paste(response_name, "~ (native_transplant + timepoint) * sym_dom + (1 | native_site:parent)")) 
  f_transplant_time_add_native_int                   <- as.formula(paste(response_name, "~ (native_transplant + timepoint) * native_site + (1 | native_site:parent)")) 
  f_sym_time_add_native_int                          <- as.formula(paste(response_name, "~ (sym_dom + timepoint) * native_site + (1 | native_site:parent)"))
  f_sym_time_add_transplant_int                     <- as.formula(paste(response_name, "~ (sym_dom + timepoint) * native_transplant + (1 | native_site:parent)"))
  # ---- 3 predictors, fully interactive, 4th predictor absent ----
  f_site_int_sym_int                                    <- as.formula(paste(response_name, "~ native_site * native_transplant * sym_dom + (1 | native_site:parent)"))
  f_site_int_time_int                                   <- as.formula(paste(response_name, "~ native_site * native_transplant * timepoint + (1 | native_site:parent)"))
  f_native_sym_time_int                                 <- as.formula(paste(response_name, "~ native_site * sym_dom * timepoint + (1 | native_site:parent)"))             
  f_transplant_sym_time_int                             <- as.formula(paste(response_name, "~ native_transplant * sym_dom * timepoint + (1 | native_site:parent)"))        
  # ---- 4 predictors, fully additive ----
  f_sym_time_site_add                                   <- as.formula(paste(response_name, "~ sym_dom + timepoint + native_site + native_transplant + (1 | native_site:parent)"))
  # ---- 4 predictors, one pair interacts + other two additive ----
  f_site_int_sym_time_add                                     <- as.formula(paste(response_name, "~ native_site * native_transplant + sym_dom + timepoint + (1 | native_site:parent)"))     
  f_site_add_sym_time_int                                     <- as.formula(paste(response_name, "~ native_site + native_transplant + timepoint * sym_dom + (1 | native_site:parent)"))
  f_native_sym_int_transplant_time_add                        <- as.formula(paste(response_name, "~ native_site * sym_dom + timepoint + native_transplant + (1 | native_site:parent)"))     
  f_transplant_sym_int_native_time_add                        <- as.formula(paste(response_name, "~ native_transplant * sym_dom + timepoint + native_site + (1 | native_site:parent)"))
  f_native_time_int_transplant_sym_add                        <- as.formula(paste(response_name, "~ native_site * timepoint + native_transplant + sym_dom + (1 | native_site:parent)"))     
  f_transplant_time_int_native_sym_add                        <- as.formula(paste(response_name, "~ native_transplant * timepoint + sym_dom + native_site + (1 | native_site:parent)"))
  # ---- 4 predictors, one triple fully interacts + 4th additive ----
  f_site_sym_int_time_add                                    <- as.formula(paste(response_name, "~ native_site * native_transplant * sym_dom + timepoint + (1 | native_site:parent)"))
  f_site_time_int_sym_add                                    <- as.formula(paste(response_name, "~ native_site * native_transplant * timepoint + sym_dom + (1 | native_site:parent)"))
  f_native_sym_time_int_transplant_add                       <- as.formula(paste(response_name, "~ native_site * sym_dom * timepoint + native_transplant + (1 | native_site:parent)"))
  f_transplant_sym_time_int_native_add                       <- as.formula(paste(response_name, "~ native_transplant * sym_dom * timepoint + native_site + (1 | native_site:parent)"))
  # ---- 4 predictors, fully interactive ----
  f_site_sym_time_int                                       <- as.formula(paste(response_name, "~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent)"))  # fixed: was missing "+"
  
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
    native_sym_int_time_add               = fit_one(f_native_sym_int_time_add, data),
    native_sym_int_transplant_add         = fit_one(f_native_sym_int_transplant_add, data),
    native_time_int_sym_add               = fit_one(f_native_time_int_sym_add, data),
    native_time_int_transplant_add        = fit_one(f_native_time_int_transplant_add, data),
    transplant_sym_int_time_add           = fit_one(f_transplant_sym_int_time_add, data),
    transplant_sym_int_native_add         = fit_one(f_transplant_sym_int_native_add, data),
    transplant_time_int_sym_add           = fit_one(f_transplant_time_int_sym_add, data),
    transplant_time_int_native_add        = fit_one(f_transplant_time_int_native_add, data),
    sym_time_int_native_add               = fit_one(f_sym_time_int_native_add, data),
    sym_time_int_transplant_add           = fit_one(f_sym_time_int_transplant_add, data),

    site_add_sym_int                      = fit_one(f_site_add_sym_int, data),                               
    site_add_time_int                     = fit_one(f_site_add_time_int, data),                         
    native_sym_add_time_int               = fit_one(f_native_sym_add_time_int, data),                        
    native_sym_add_transplant_int         = fit_one(f_native_sym_add_transplant_int, data),                  
    native_time_add_sym_int               = fit_one(f_native_time_add_sym_int, data),            
    native_time_add_transplant_int        = fit_one(f_native_time_add_transplant_int, data),                  
    transplant_sym_add_time_int           = fit_one(f_transplant_sym_add_time_int, data),           
    transplant_sym_add_native_int         = fit_one(f_transplant_sym_add_native_int, data),                       
    transplant_time_add_sym_int           = fit_one(f_transplant_time_add_sym_int, data),            
    transplant_time_add_native_int        = fit_one(f_transplant_time_add_native_int, data),              
    sym_time_add_native_int               = fit_one(f_sym_time_add_native_int, data),           
    sym_time_add_transplant_int           = fit_one(f_sym_time_add_transplant_int, data),                 

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

#fit models within each timepoint
fit_models_within_timepoint <- function(response_name, data) {
  # ---- 0 predictors ----
  f_null                 <- as.formula(paste(response_name, "~ 1 + (1 | native_site:parent)"))
  # ---- 1 predictor ----
  f_native                <- as.formula(paste(response_name, "~ native_site + (1 | native_site:parent)"))
  f_transplant            <- as.formula(paste(response_name, "~ native_transplant + (1 | native_site:parent)"))
  f_sym                   <- as.formula(paste(response_name, "~ sym_dom + (1 | native_site:parent)"))
  # ---- 2 predictors, additive ----
  f_site_add                  <- as.formula(paste(response_name, "~ native_site + native_transplant + (1 | native_site:parent)"))
  f_native_sym                <- as.formula(paste(response_name, "~ native_site + sym_dom + (1 | native_site:parent)"))
  f_transplant_sym            <- as.formula(paste(response_name, "~ native_transplant + sym_dom + (1 | native_site:parent)"))
  # ---- 2 predictors, interactive ----
  f_site_int                  <- as.formula(paste(response_name, "~ native_site * native_transplant + (1 | native_site:parent)"))
  f_native_sym_int            <- as.formula(paste(response_name, "~ native_site * sym_dom + (1 | native_site:parent)"))
  f_transplant_sym_int        <- as.formula(paste(response_name, "~ native_transplant * sym_dom + (1 | native_site:parent)"))
  # ---- 3 predictors, all additive ----
  f_site_add_sym_add          <- as.formula(paste(response_name, "~ native_site + native_transplant + sym_dom + (1 | native_site:parent)"))
  # ---- 3 predictors, one pair interacts + third additive ----
  f_site_int_sym_add                          <- as.formula(paste(response_name, "~ native_site * native_transplant + sym_dom + (1 | native_site:parent)"))
  f_native_sym_int_transplant_add             <- as.formula(paste(response_name, "~ native_site * sym_dom + native_transplant + (1 | native_site:parent)"))
  f_transplant_sym_int_native_add             <- as.formula(paste(response_name, "~ native_transplant * sym_dom + native_site + (1 | native_site:parent)"))
  # ---- 3 predictors, "hub" (one variable interacts with each of the
  # other two separately, those two don't interact with each other) ----
  f_site_add_sym_int                                 <- as.formula(paste(response_name, "~ (native_site + native_transplant) * sym_dom + (1 | native_site:parent)"))
  f_native_sym_add_transplant_int                    <- as.formula(paste(response_name, "~ (native_site + sym_dom) * native_transplant + (1 | native_site:parent)"))
  f_transplant_sym_add_native_int                    <- as.formula(paste(response_name, "~ (native_transplant + sym_dom) * native_site + (1 | native_site:parent)"))          
  # ---- 3 predictors, fully interactive ----
  f_site_int_sym_int                                    <- as.formula(paste(response_name, "~ native_site * native_transplant * sym_dom + (1 | native_site:parent)"))

  models <- list(
    null                                  = fit_one(f_null, data),
    
    native                                = fit_one(f_native, data),
    transplant                            = fit_one(f_transplant, data),
    sym                                   = fit_one(f_sym, data),
    
    site_add                              = fit_one(f_site_add, data),
    native_sym                            = fit_one(f_native_sym, data),
    transplant_sym                        = fit_one(f_transplant_sym, data),

    site_int                              = fit_one(f_site_int, data),
    native_sym_int                        = fit_one(f_native_sym_int, data),
    transplant_sym_int                    = fit_one(f_transplant_sym_int, data),

    site_add_sym_add                      = fit_one(f_site_add_sym_add, data),
    
    site_int_sym_add                      = fit_one(f_site_int_sym_add, data),                         
    native_sym_int_transplant_add         = fit_one(f_native_sym_int_transplant_add, data),
    transplant_sym_int_native_add         = fit_one(f_transplant_sym_int_native_add, data),
    
    site_add_sym_int                      = fit_one(f_site_add_sym_int, data),                               
    native_sym_add_transplant_int         = fit_one(f_native_sym_add_transplant_int, data),                  
    transplant_sym_add_native_int         = fit_one(f_transplant_sym_add_native_int, data),                       

    site_int_sym_int                      = fit_one(f_site_int_sym_int, data)
  )
  
  return(models)
}

#fit models at T0
fit_models_T0 <- function(response_name, data) {
  # ---- 0 predictors ----
  f_null                 <- as.formula(paste(response_name, "~ 1 + (1 | native_site:parent)"))
  # ---- 1 predictor ----
  f_native                <- as.formula(paste(response_name, "~ native_site + (1 | native_site:parent)"))
  f_sym                   <- as.formula(paste(response_name, "~ sym_dom + (1 | native_site:parent)"))
  # ---- 2 predictors, additive ----
  f_native_sym                <- as.formula(paste(response_name, "~ native_site + sym_dom + (1 | native_site:parent)"))
  # ---- 2 predictors, interactive ----
  f_native_sym_int            <- as.formula(paste(response_name, "~ native_site * sym_dom + (1 | native_site:parent)"))
   
  models <- list(
    null                                  = fit_one(f_null, data),
    
    native                                = fit_one(f_native, data),
    sym                                   = fit_one(f_sym, data),
    
    native_sym                            = fit_one(f_native_sym, data),

    native_sym_int                        = fit_one(f_native_sym_int, data)
  )
  
  return(models)
}


##compare models
compare_models_by_aic <- function(models_list, response_name) {
  
  good_models <- models_list[!sapply(models_list, is.null)]
  
  # If EVERY candidate model failed, return an empty table instead of erroring.
  if (length(good_models) == 0) {
    return(data.frame(
      response = character(0), model = character(0),
      df = numeric(0), AIC = numeric(0), singular = logical(0),
      deltaAIC = numeric(0), weight = numeric(0)
    ))
  }
  
  model_names     <- names(good_models)
  aic_values      <- sapply(good_models, AIC)
  df_values       <- sapply(good_models, function(m) attr(logLik(m), "df"))
  singular_values <- sapply(good_models, isSingular)
  
  results <- data.frame(
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

dir.create(here("Output", "Physiology"), recursive = TRUE, showWarnings = FALSE)

##### Siderastrea siderea models across time #####
sym_physio_sid <- sym_physio %>% filter(coral_species == "Siderastrea siderea")

#variables to use for sid models
#"chla_log", "chla_sym_log", "sym_density_log", "dry_weight", "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs_log", "calcification" 

response_vars_sid <- c("chla_log", "chla_sym_log", "sym_density_log", "dry_weight", 
                             "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs_log", "calcification")

#create objects to put stuff in for the models
all_fitted_models_sid <- list()
all_aic_tables_sid    <- list()
best_models_summary_sid <- tibble(
  response = character(),
  best_model = character(),
  model_formula = character(),
  AIC = double(),
  significant_predictors = character()
)

#run models for all of these variables adding in symbiont type to see what the best predictor is
#siderastrea siderea
for (response_name in response_vars_sid) {
  
  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models_sid[[response_name]] <- list()
  all_aic_tables_sid[[response_name]]    <- list()
  
  this_data <- sym_physio_sid %>% 
    filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
    droplevels()
    
    cat("Rows of data after filtering:", nrow(this_data), "\n")
    
    if (nrow(this_data) < 10) {
      cat("Not enough data - skipping.\n")
      next
    }
    
    # ---- FIT all candidate models ----
    models_list <- fit_models_with_time(response_name, this_data)
    all_fitted_models_sid[[response_name]] <- models_list
    
    n_fit    <- sum(!sapply(models_list, is.null))
    n_failed <- sum(sapply(models_list, is.null))
    cat(n_fit, "of models fit successfully;", n_failed, "failed.\n")
    
    # ---- AIC COMPARE ----
    aic_table <- compare_models_by_aic(models_list, response_name)
    all_aic_tables_sid[[response_name]]<- aic_table
    
    if (nrow(aic_table) == 0) {
      cat("None of the models fit for", response_name, "- see [ERROR] messages above.\n")
      next
    }
    
    print(aic_table)
    
    file_tag <- paste0(response_name)
    write_csv(aic_table, here("Output", "Physiology", paste0("aic_table_", file_tag, "_sid.csv")))
    
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
        pdf(here("Output", "Physiology", paste0("checkmodel_", file_tag, "_sid.pdf")), width = 12, height = 10)
        print(check_model(best_model_reml))
      }, error = function(e) {
        cat("check_model() failed:", conditionMessage(e), "\n")
      }, finally = {
        dev.off()
      })
    }
    # ---- Pull significant predictors (p < 0.05) from the best model ----
    significant_terms <- NA_character_
    if (!is.null(best_model_reml)) {
      anova_result <- tryCatch(car::Anova(best_model_reml), error = function(e) NULL)
      if (!is.null(anova_result)) {
        sig_rows <- rownames(anova_result)[anova_result[["Pr(>Chisq)"]] < 0.05]
        significant_terms <- if (length(sig_rows) > 0) paste(sig_rows, collapse = ", ") else "none"
      }
    }
    
    # ---- Record this trait's result in the summary table ----
    best_models_summary_sid <- best_models_summary_sid %>%
      add_row(
        response = response_name,
        best_model = best_model_name,
        model_formula = if (!is.null(best_model_reml)) {
          paste(deparse(formula(best_model_reml)), collapse = " ")
        } else {
          NA_character_
        },
        AIC = aic_table$AIC[1],
        significant_predictors = significant_terms
      )
}

#look at best models summary for siderastrea siderea over time
#print(best_models_summary_sid)
write_csv(best_models_summary_sid, here("Output", "Physiology", "best_models_summary_time_sid.csv"))

#### Siderastrea siderea models at each timepoint separately #####
#separate data into timepoints
#just switched out T4 and T12 for this section
sym_physio_sid_T4 <- sym_physio_sid %>% filter(timepoint == "T4")

#create objects to put stuff in for the models
all_fitted_models_sid_T4 <- list()
all_aic_tables_sid_T4    <- list()
best_models_summary_sid_T4 <- tibble(
  response = character(),
  best_model = character(),
  model_formula = character(),
  AIC = double(),
  significant_predictors = character()
)

#run models for all of these variables adding in symbiont type to see what the best predictor is
#siderastrea siderea
for (response_name in response_vars_sid) {
  
  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models_sid_T4[[response_name]] <- list()
  all_aic_tables_sid_T4[[response_name]]    <- list()
  
  this_data <- sym_physio_sid_T4 %>% 
    filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
    droplevels()
  
  cat("Rows of data after filtering:", nrow(this_data), "\n")
  
  if (nrow(this_data) < 10) {
    cat("Not enough data - skipping.\n")
    next
  }
  
  # ---- FIT all candidate models ----
  models_list <- fit_models_within_timepoint(response_name, this_data)
  all_fitted_models_sid_T4[[response_name]] <- models_list
  
  n_fit    <- sum(!sapply(models_list, is.null))
  n_failed <- sum(sapply(models_list, is.null))
  cat(n_fit, "of models fit successfully;", n_failed, "failed.\n")
  
  # ---- AIC COMPARE ----
  aic_table <- compare_models_by_aic(models_list, response_name)
  all_aic_tables_sid_T4[[response_name]]<- aic_table
  
  if (nrow(aic_table) == 0) {
    cat("None of the models fit for", response_name, "- see [ERROR] messages above.\n")
    next
  }
  
  print(aic_table)
  
  file_tag <- paste0(response_name)
  write_csv(aic_table, here("Output", "Physiology", paste0("aic_table_", file_tag, "_sid_T4.csv")))
  
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
      pdf(here("Output", "Physiology", paste0("checkmodel_", file_tag, "_sid_T4.pdf")), width = 12, height = 10)
      print(check_model(best_model_reml))
    }, error = function(e) {
      cat("check_model() failed:", conditionMessage(e), "\n")
    }, finally = {
      dev.off()
    })
  }
  # ---- Pull significant predictors (p < 0.05) from the best model ----
  significant_terms <- NA_character_
  if (!is.null(best_model_reml)) {
    anova_result <- tryCatch(car::Anova(best_model_reml), error = function(e) NULL)
    if (!is.null(anova_result)) {
      sig_rows <- rownames(anova_result)[anova_result[["Pr(>Chisq)"]] < 0.05]
      significant_terms <- if (length(sig_rows) > 0) paste(sig_rows, collapse = ", ") else "none"
    }
  }
  
  # ---- Record this trait's result in the summary table ----
  best_models_summary_sid_T4 <- best_models_summary_sid_T4 %>%
    add_row(
      response = response_name,
      best_model = best_model_name,
      model_formula = if (!is.null(best_model_reml)) {
        paste(deparse(formula(best_model_reml)), collapse = " ")
      } else {
        NA_character_
      },
      AIC = aic_table$AIC[1],
      significant_predictors = significant_terms
    )
}

#look at best models summary for siderastrea siderea over time
#print(best_models_summary_sid_T4)
write_csv(best_models_summary_sid_T4, here("Output", "Physiology", "best_models_summary_time_sid_T4.csv"))

##### Branching Porites sp. models across time ####
sym_physio_por <- sym_physio %>% filter(coral_species == "Porites sp.")

#variables to use for por models
#"chla_sqrt", "chla_sym_sqrt", "sym_density_sqrt", "dry_weight", "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs", "calcification" 

response_vars_por <- c("chla_sqrt", "chla_sym_sqrt", "sym_density_sqrt", "dry_weight", 
                       "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs", "calcification" )

#create objects to put stuff in for the models
all_fitted_models_por <- list()
all_aic_tables_por    <- list()
best_models_summary_por <- tibble(
  response = character(),
  best_model = character(),
  model_formula = character(),
  AIC = double(),
  significant_predictors = character()
)

#run models for all of these variables adding in symbiont type to see what the best predictor is
#Branching Porites sp.
for (response_name in response_vars_por) {
  
  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models_por[[response_name]] <- list()
  all_aic_tables_por[[response_name]]    <- list()
  
  this_data <- sym_physio_por %>% 
    filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
    droplevels()
  
  cat("Rows of data after filtering:", nrow(this_data), "\n")
  
  if (nrow(this_data) < 10) {
    cat("Not enough data - skipping.\n")
    next
  }
  
  # ---- FIT all candidate models ----
  models_list <- fit_models_with_time(response_name, this_data)
  all_fitted_models_por[[response_name]] <- models_list
  
  n_fit    <- sum(!sapply(models_list, is.null))
  n_failed <- sum(sapply(models_list, is.null))
  cat(n_fit, "of models fit successfully;", n_failed, "failed.\n")
  
  # ---- AIC COMPARE ----
  aic_table <- compare_models_by_aic(models_list, response_name)
  all_aic_tables_por[[response_name]]<- aic_table
  
  if (nrow(aic_table) == 0) {
    cat("None of the models fit for", response_name, "- see [ERROR] messages above.\n")
    next
  }
  
  print(aic_table)
  
  file_tag <- paste0(response_name)
  write_csv(aic_table, here("Output", "Physiology", paste0("aic_table_", file_tag, "_por.csv")))
  
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
      pdf(here("Output", "Physiology", paste0("checkmodel_", file_tag, "_por.pdf")), width = 12, height = 10)
      print(check_model(best_model_reml))
    }, error = function(e) {
      cat("check_model() failed:", conditionMessage(e), "\n")
    }, finally = {
      dev.off()
    })
  }
  # ---- Pull significant predictors (p < 0.05) from the best model ----
  significant_terms <- NA_character_
  if (!is.null(best_model_reml)) {
    anova_result <- tryCatch(car::Anova(best_model_reml), error = function(e) NULL)
    if (!is.null(anova_result)) {
      sig_rows <- rownames(anova_result)[anova_result[["Pr(>Chisq)"]] < 0.05]
      significant_terms <- if (length(sig_rows) > 0) paste(sig_rows, collapse = ", ") else "none"
    }
  }
  
  # ---- Record this trait's result in the summary table ----
  best_models_summary_por <- best_models_summary_por %>%
    add_row(
      response = response_name,
      best_model = best_model_name,
      model_formula = if (!is.null(best_model_reml)) {
        paste(deparse(formula(best_model_reml)), collapse = " ")
      } else {
        NA_character_
      },
      AIC = aic_table$AIC[1],
      significant_predictors = significant_terms
    )
}

#look at best models summary for Branching Porites sp. over time
#print(best_models_summary_por)
write_csv(best_models_summary_por, here("Output", "Physiology", "best_models_summary_time_por.csv"))

#### Branching Porites sp. models at each timepoint separately #####
#separate data into timepoints
#just switched out T4 and T12 for this section
sym_physio_por_T4 <- sym_physio_por %>% filter(timepoint == "T4")

#create objects to put stuff in for the models
all_fitted_models_por_T4 <- list()
all_aic_tables_por_T4    <- list()
best_models_summary_por_T4 <- tibble(
  response = character(),
  best_model = character(),
  model_formula = character(),
  AIC = double(),
  significant_predictors = character()
)

#run models for all of these variables adding in symbiont type to see what the best predictor is
#Branching Porites sp.
for (response_name in response_vars_por) {
  
  cat("\n\n=====", response_name, "=====\n")
  all_fitted_models_por_T4[[response_name]] <- list()
  all_aic_tables_por_T4[[response_name]]    <- list()
  
  this_data <- sym_physio_por_T4 %>% 
    filter(!is.na(.data[[response_name]]), !is.na(sym_dom)) %>%
    droplevels()
  
  cat("Rows of data after filtering:", nrow(this_data), "\n")
  
  if (nrow(this_data) < 10) {
    cat("Not enough data - skipping.\n")
    next
  }
  
  # ---- FIT all candidate models ----
  models_list <- fit_models_within_timepoint(response_name, this_data)
  all_fitted_models_por_T4[[response_name]] <- models_list
  
  n_fit    <- sum(!sapply(models_list, is.null))
  n_failed <- sum(sapply(models_list, is.null))
  cat(n_fit, "of models fit successfully;", n_failed, "failed.\n")
  
  # ---- AIC COMPARE ----
  aic_table <- compare_models_by_aic(models_list, response_name)
  all_aic_tables_por_T4[[response_name]]<- aic_table
  
  if (nrow(aic_table) == 0) {
    cat("None of the models fit for", response_name, "- see [ERROR] messages above.\n")
    next
  }
  
  print(aic_table)
  
  file_tag <- paste0(response_name)
  write_csv(aic_table, here("Output", "Physiology", paste0("aic_table_", file_tag, "_por_T4.csv")))
  
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
      pdf(here("Output", "Physiology", paste0("checkmodel_", file_tag, "_por_T4.pdf")), width = 12, height = 10)
      print(check_model(best_model_reml))
    }, error = function(e) {
      cat("check_model() failed:", conditionMessage(e), "\n")
    }, finally = {
      dev.off()
    })
  }
  # ---- Pull significant predictors (p < 0.05) from the best model ----
  significant_terms <- NA_character_
  if (!is.null(best_model_reml)) {
    anova_result <- tryCatch(car::Anova(best_model_reml), error = function(e) NULL)
    if (!is.null(anova_result)) {
      sig_rows <- rownames(anova_result)[anova_result[["Pr(>Chisq)"]] < 0.05]
      significant_terms <- if (length(sig_rows) > 0) paste(sig_rows, collapse = ", ") else "none"
    }
  }
  
  # ---- Record this trait's result in the summary table ----
  best_models_summary_por_T4 <- best_models_summary_por_T4 %>%
    add_row(
      response = response_name,
      best_model = best_model_name,
      model_formula = if (!is.null(best_model_reml)) {
        paste(deparse(formula(best_model_reml)), collapse = " ")
      } else {
        NA_character_
      },
      AIC = aic_table$AIC[1],
      significant_predictors = significant_terms
    )
}

#look at best models summary for branching porites sp. over time
#print(best_models_summary_por_T4)
write_csv(best_models_summary_por_T4, here("Output", "Physiology", "best_models_summary_time_por_T4.csv"))


##### Simple plots of all variables to look at all data #####

##### Exploratory plots #####
#chlc2	chla	chlc2_sym	chla_sym	dry_weight	sym_density	calcification	pct_change_dw	pct_change_sa 
#create a function to plot all!!

#species
ss_data <- sym_physio %>% filter(coral_species == "Siderastrea siderea") %>% droplevels()
pp_data <- sym_physio %>% filter(coral_species == "Porites sp.") %>% droplevels()
species_data <- list("Siderastrea siderea" = ss_data,"Porites sp." = pp_data)

#variables to examine
response_vars_plots <- c("chla", "chla_sym", "dry_weight", "sym_density","calcification", "pct_change_dw", "pct_change_sa",
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

#Physio site time sym
iwalk(species_data, function(d, sp) {
  walk(response_vars_plots, function(resp) {
    p <- plot_physio_site_time_sym(d, resp)
    ggsave(
      here("Output", "Physiology", "Plots", paste0(sp, "_", resp, "site_time_sym.pdf")),
      p, height = 5, width = 12
    )
  })
})

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
iwalk(species_data, function(d, sp) {
  walk(response_vars_plots, function(resp) {
    p <- plot_physio_sym(d, resp)
    ggsave(
      here("Output", "Physiology", "Plots", paste0(sp, "_", resp, "_sym.pdf")),
      p, height = 5, width = 8
    )
  })
})


##### Final best fit models and pairwise comparisons to add to plots #####

###Siderastrea siderea
#"chla_log", "chla_sym_log", "sym_density_log", "dry_weight", 
#"PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs_log", "calcification"

###all models across time: models with a significant effect of symbiont identity to examine further
##only affects bulk physiology
#sid_mod_PR <- lmer(PGROSS.RDAY	~ native_transplant * timepoint + native_site + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#sid_mod_PNET <- lmer(PNET	~ native_transplant * timepoint + native_site + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#sid_mod_PGROSS <- lmer(PGROSS	~ native_transplant * timepoint + native_site + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#sid_mod_RDARK <- lmer(RDARK_abs	~ native_site + timepoint + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#sid_mod_LEDR <- lmer(LEDR_abs_log	~ (native_site + native_transplant) * timepoint + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_chla <- lmer(chla_log	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_chla <- lmer(chla_log ~ sym_dom + (1|native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
summary(sid_mod_chla)
anova(sid_mod_chla)
#significant terms: 
#native_transplant, timepoint, native_site:timepoint, native_site:sym_dom, native_transplant:timepoint, timepoint:sym_dom
ranova(sid_mod_chla) #ns

sid_emm_chla <-emmeans(sid_mod_chla, pairwise ~  sym_dom) 
sid_emm_chla <- emmeans(sid_mod_chla, specs = c("sym_dom"))
pairs(sid_emm_chla)

sid_mod_chla_sym <- lmer(chla_sym_log	~ native_site + sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_sym <- lmer(sym_density_log	~ native_site * sym_dom * timepoint + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_dw <- lmer(dry_weight	~ native_site * native_transplant * timepoint * sym_dom + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_calc <- lmer(calcification	~ native_site * sym_dom + timepoint + (1 | native_site:parent), data = sym_physio_sid, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))

##models at T12 with symbiont identity included to examine further
#same as general models but also including PNET
sid_mod_T12_chla <- lmer(chla_log	~ (native_site + native_transplant) * sym_dom + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_T12_chla_sym <- lmer(chla_sym_log	~ native_site * sym_dom + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_T12_sym <- lmer(sym_density_log	~ native_site + sym_dom + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_T12_dw <- lmer(dry_weight	~ ~ (native_transplant + sym_dom) * native_site + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_T12_calc <- lmer(calcification	~ native_site * sym_dom + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
sid_mod_T12_PNET <- lmer(PNET	~ native_site + sym_dom + (1 | native_site:parent), data = sym_physio_sid_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))

### Branching Porites sp.
#"chla_sqrt", "chla_sym_sqrt", "sym_density_sqrt", "dry_weight", "PGROSS.RDAY", "PNET", "PGROSS", "RDARK_abs", "LEDR_abs", "calcification"
#normality of all models checked using the check_model from above which prints for best model
#can be found in Output, Physiology

###all models across time: models with a significant effect of symbiont identity
#por_mod_chla <- lmer(chla_sqrt	~ native_site * native_transplant * timepoint + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#por_mod_chla_sym <- lmer(chla_sym_sqrt	~ native_transplant + timepoint + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#por_mod_sym <- lmer(sym_density_sqrt	~ native_site * native_transplant + sym_dom + timepoint + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#por_mod_dw <- lmer(dry_weight	~ native_site * native_transplant * timepoint + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#por_mod_PR <- lmer(PGROSS.RDAY	~ native_transplant + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#por_mod_PGROSS <- lmer(PGROSS	~ native_site * timepoint + native_transplant + sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
#RDARK_abs	NULL
#por_mod_calc <- lmer(calcification	~ (native_site + timepoint) * native_transplant + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
por_mod_PNET <- lmer(PNET	~ native_site * timepoint + native_transplant + sym_dom + (1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
por_mod_LEDR <- lmer(LEDR_abs	~ native_transplant * timepoint + sym_dom +(1 | native_site:parent), data = sym_physio_por, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))

##models at T12 with symbiont identity included to examine further
por_mod_T12_RDARK <- lmer(RDARK_abs	~ native_site + sym_dom +(1 | native_site:parent), data = sym_physio_por_T12, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))

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


#### Ordination plots of physiology data ####
#MP just doing this to look at it to see if anything interesting pops out

#ordination of physio data
#load physio data and generate dataframe with id as rownames
#dataframe = sym_physio
sym_physio_data <- sym_physio %>% select(sample_id, chla, chla_sym, dry_weight, sym_density, PNET, PGROSS, RDARK, LEDR, PGROSS.RDAY) %>% 
  column_to_rownames("sample_id") %>% drop_na()
sym_physio_data_T12 <- sym_physio %>% filter(timepoint == "T12") %>%
  select(sample_id, chla, chla_sym, dry_weight, sym_density, calcification, PNET, PGROSS, RDARK, LEDR, PGROSS.RDAY) %>% 
  column_to_rownames("sample_id") %>% drop_na()

#make them matrices
sym_physio_mat <- as.matrix(sym_physio_data)
sym_physio_mat_T12 <- as.matrix(sym_physio_data_T12)
#quick glance at data to look for strong patterns
pairs(x = sym_physio_mat, gap = 0, cex.labels = 0.5) #look similar even between metrics
pairs(x = sym_physio_mat_T12, gap = 0, cex.labels = 0.5)
#scale across variable types
sym_physio_mat <- scale(sym_physio_mat)
sym_physio_mat_T12 <- scale(sym_physio_mat_T12)
#generate pca data
pca_sym_physio <- prcomp(sym_physio_mat)
pca_sym_physio_T12 <- prcomp(sym_physio_mat_T12)
#collapse PCA data
pc_axes_physio <- as.data.frame(pca_sym_physio$x)
pc_axes_physio_T12 <- as.data.frame(pca_sym_physio_T12$x)
#add sample_id back
pc_axes_physio$sample_id <- rownames(pc_axes_physio) 
pc_axes_physio_T12$sample_id <- rownames(pc_axes_physio_T12) 

#create metadata
phys_meta <- sym_physio %>% select(sample_id, coral_species, transplant_group, native_site, 
                                   native_transplant, transplant_group_sampled, sym_dom)
#put data back with metadata
pca_physio <- pc_axes_physio %>% left_join(phys_meta, by = "sample_id")
pca_physio_T12 <- pc_axes_physio_T12 %>% left_join(phys_meta, by = "sample_id")

#plot sym across all timepoints
pca_physio_sym <- ggplot(pca_physio, aes(x = PC1, y = PC2, color = sym_dom)) +
  geom_point(size = 3) +
  scale_color_manual(values = sym_dom_colors) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)
  stat_ellipse(aes(group = sym_dom), level = 0.95, alpha = 0.5, linewidth = 0.8) +
  theme_classic(base_size = 22) +
  facet_wrap(coral_species~., scales = "free")
pca_physio_sym
#plot site across all timepoints
pca_physio_site <- ggplot(pca_physio, aes(x = PC1, y = PC2, color = transplant_group)) +
  geom_point(size = 3) +
  scale_color_manual(values = site_colors) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)
  stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8) +
  theme_classic(base_size = 22) +
  facet_wrap(coral_species ~., scales = "free")
pca_physio_site

#plot sym at T12
pca_physio_sym_T12 <- ggplot(pca_physio_T12, aes(x = PC1, y = PC2, color = sym_dom)) +
  geom_point(size = 3) +
  scale_color_manual(values = sym_dom_colors) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)
  stat_ellipse(aes(group = sym_dom), level = 0.95, alpha = 0.5, linewidth = 0.8) +
  theme_classic(base_size = 22) +
  facet_wrap(coral_species~., scales = "free")
pca_physio_sym_T12
#plot site at T12
pca_physio_site_T12 <- ggplot(pca_physio_T12, aes(x = PC1, y = PC2, color = transplant_group)) +
  geom_point(size = 3) +
  scale_color_manual(values = site_colors) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)
  stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8) +
  theme_classic(base_size = 22) +
  facet_wrap(coral_species ~., scales = "free")
pca_physio_site_T12

#don't really see much of interest here - maybe more interesting if we separate bulk physio from respirometry??

#ggsave(here("Output", "Physiology", "physio_params_T12_pc_spp_perf.pdf"), pca_physio_spp, h = 5, w = 8)

# physio_fit <- envfit(pca_sym_physio$x[, c("PC1", "PC2")], sym_physio_data, permutations = 999, na.rm = TRUE)
# physio_T12_fit <- envfit(pca_sym_physio_T12$x[, c("PC1", "PC2")], sym_physio_data_T12, permutations = 999, na.rm = TRUE)
# 
# physio_scores <- as.data.frame(scores(physio_fit, display = "vectors")) %>% mutate(variable = rownames(.))
# physio_T12_scores <- as.data.frame(scores(physio_T12_fit, display = "vectors")) %>% mutate(variable = rownames(.))
# 
# pca_physio_spp_arrows <- ggplot(pca_physio, aes(x = PC1, y = PC2, color = morphoT12y, fill = morphoT12y)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = morph_colors, name = "MorphoT12y", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
#   scale_fill_manual(values = morph_colors, name = "MorphoT12y", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
#   #scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
#   #scale_fill_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
#   stat_ellipse(geom = "polygon", alpha = 0.1)+
#   #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)+
#   #stat_ellipse(aes(fill = group), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
#   geom_segment(data = physio_scores,
#                aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
#                arrow = arrow(length = unit(0.25, "cm")),
#                color = "black", inherit.aes = FALSE) +
#   #geom_text(data = physio_scores,
#   #          aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
#   #          color = "black", size = 5, inherit.aes = FALSE) +
#   theme_classic(base_size = 22)
# 
# pca_physio_spp_arrows
# 
# ggsave(here("Output", "PhysioT12y", "physio_params_morph_arrows_nolabels.pdf"), pca_physio_spp_arrows, h = 8, w = 12)

