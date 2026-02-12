library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)

# Source your fitting functions
source("functions/Reco_Temp_GPP_PAR_fitting_functions.R")

selected_year <- 2024
# ============================================================================
# 1. LOAD AND PREPARE FLUXNET DATA
# ============================================================================

# Load FLUXNET data
fluxnet_data <- fread("C:/github/GHG-data-useful-scripts/data/ICOSETC_SE-Deg_FLUXNET_INTERIM_HH_L2/ICOSETC_SE-Deg_FLUXNET_INTERIM_HH_L2.csv")

# Prepare data

ec_data_sub <- fluxnet_data %>%
  # Replace -9999 with NA across entire dataset
  mutate(across(everything(), ~na_if(., -9999))) %>%
  mutate(
    datetime = ymd_hm(as.character(TIMESTAMP_END)),
    Temp = TS_F_MDS_1, # or TA_F for air temp
    PARin = PPFD_IN,
    Reco = RECO_NT_VUT_REF,
    GPP = GPP_NT_VUT_REF
  ) %>%
  select(datetime, Temp, PARin, Reco, GPP) %>%
  filter(year(datetime)==selected_year) #all year
# filter(year(datetime)==selected_year, PARin < 0.5) # All year, only nighttime

# ============================================================================
# 2. FIT RECO TEMPERATURE RESPONSE
# ============================================================================

cat("\n=== Fitting Reco temperature response ===\n")

# Fit Lloyd & Taylor to existing Reco data
reco_fit <- fit_reco_lloyd_taylor(
  temp = ec_data_sub$Temp,
  reco = ec_data_sub$Reco
)

# Plot
plot_reco_response(reco_fit)

# ============================================================================
# 3. FIT GPP LIGHT RESPONSE
# ============================================================================

ec_data_sub <- fluxnet_data %>%
  # Replace -9999 with NA across entire dataset
  mutate(across(everything(), ~na_if(., -9999))) %>%
  mutate(
    datetime = ymd_hm(as.character(TIMESTAMP_END)),
    Temp = TS_F_MDS_3, # or TA_F for air temp
    PARin = PPFD_IN,
    Reco = RECO_NT_VUT_REF,
    GPP = GPP_NT_VUT_REF
  ) %>%
  select(datetime, Temp, PARin, Reco, GPP) %>%
  filter(year(datetime)==selected_year, month(datetime) %in% 6:8) # Only peak summer for the GPP curves maybe

cat("\n=== Fitting GPP light response ===\n")

# Fit hyperbolic light response to existing GPP data
# Note: FLUXNET GPP is positive, so negate it
GPPit <- fit_gpp_hyperbolic(
  par = ec_data_sub$PARin,
  gpp = -ec_data_sub$GPP,  # (provide negative GPP values)
  par_min_threshold=10
)

# Plot
plot_gpp_response(GPPit)

# ============================================================================
# 4. SAVE RESULTS
# ============================================================================

# Create summary dataframe [You can use something like this and edit a bit the codes to save all years separately]
results_summary <- data.frame(
  Site = "SE-Deg",
  R10 = reco_fit$R10,
  R10_se = reco_fit$R10_se,
  E0 = reco_fit$E0,
  E0_se = reco_fit$E0_se,
  Reco_R2 = reco_fit$R2,
  Reco_n = reco_fit$n_obs,
  alpha = GPPit$alpha,
  Pmax = GPPit$Pmax,
  GPP_R2 = GPPit$R2,
  GPP_n = GPPit$n_obs
)

print(results_summary)

cat("\n=== Analysis complete! ===\n")
