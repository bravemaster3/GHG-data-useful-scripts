library(ggplot2)

# Step 1: Fit Reco from nighttime NEE
fit_reco_from_nee <- function(temp, nee, par, par_threshold = 0.5) {
  # temp: temperature (°C)
  # nee: net ecosystem exchange (µmol CO2 m-2 s-1)
  # par: photosynthetically active radiation (µmol m-2 s-1)
  # par_threshold: PAR threshold for nighttime (default 0.5)

  # MATLAB filtering: nighttime only (PAR < threshold), positive NEE
  # At night, NEE ≈ Reco (no photosynthesis)
  valid <- !is.na(temp) & !is.na(nee) & !is.na(par) &
    nee > 0 & par < par_threshold

  temp_clean <- temp[valid]
  nee_night <- nee[valid]

  cat("Fitting Reco from", length(temp_clean), "nighttime NEE points\n")

  # MATLAB parameterization
  fit <- nls(nee_night ~ r_ref * exp(e0 * (1/56.02 - 1/(temp_clean + 46.02))),
             start = list(r_ref = 1.2, e0 = 337),
             control = nls.control(maxiter = 1000, warnOnly = TRUE))

  params <- coef(fit)
  conf_int <- tryCatch(confint(fit), error = function(e) NULL)

  # Predictions
  temp_seq <- seq(min(temp_clean), max(temp_clean), length.out = 100)
  reco_pred <- params["r_ref"] * exp(params["e0"] *
                                       (1/56.02 - 1/(temp_seq + 46.02)))

  r_squared <- 1 - sum(residuals(fit)^2) / sum((nee_night - mean(nee_night))^2)

  result <- list(
    R_ref = params["r_ref"],
    E0 = params["e0"],
    R2 = r_squared,
    temp_data = temp_clean,
    reco_data = nee_night,
    temp_pred = temp_seq,
    reco_pred = as.vector(reco_pred),
    n_obs = length(temp_clean),
    model = fit  # Store model for predicting Reco
  )

  if (!is.null(conf_int)) {
    result$R_ref_se <- (conf_int["r_ref", 2] - conf_int["r_ref", 1]) / (2 * 1.96)
    result$E0_se <- (conf_int["e0", 2] - conf_int["e0", 1]) / (2 * 1.96)
  }

  return(result)
}


# Step 2: Estimate GPP from NEE using fitted Reco model
# GPP = Reco - NEE
estimate_gpp_from_nee <- function(temp, nee, reco_fit) {
  # temp: temperature for all time points (°C)
  # nee: NEE for all time points (µmol CO2 m-2 s-1)
  # reco_fit: result from fit_reco_from_nee()

  # Predict Reco for all temperatures using the fitted model
  params <- coef(reco_fit$model)
  reco_est <- params["r_ref"] * exp(params["e0"] *
                                      (1/56.02 - 1/(temp + 46.02)))

  # Calculate GPP: GPP = Reco - NEE
  gpp_est <- reco_est - nee

  return(data.frame(
    temp = temp,
    nee = nee,
    reco_est = reco_est,
    gpp_est = gpp_est
  ))
}


# Step 3: Fit GPP light response curve
# Corrected: Fit GPP light response using ALL estimated GPP (no PAR filtering)
fit_gpp_lightresponse <- function(par, gpp) {
  # par: photosynthetically active radiation (µmol m-2 s-1)
  # gpp: estimated GPP for ALL time points (µmol CO2 m-2 s-1)
  # NO par_threshold - uses all valid data like MATLAB

  # MATLAB filtering: only remove NaN
  valid <- !is.na(par) & !is.na(gpp)

  par_clean <- par[valid]
  gpp_clean <- gpp[valid]

  cat("Fitting GPP light response from", length(par_clean), "data points (all times)\n")

  # MATLAB parameterization: Michaelis-Menten
  fit <- nls(gpp_clean ~ alpha * par_clean * pmax / (alpha * par_clean + pmax),
             start = list(alpha = -0.1, pmax = -25),
             control = nls.control(maxiter = 1000, warnOnly = TRUE))

  params <- coef(fit)
  conf_int <- tryCatch(confint(fit), error = function(e) NULL)

  # Predictions
  par_seq <- seq(0, max(par_clean), length.out = 100)
  gpp_pred <- params["alpha"] * par_seq * params["pmax"] /
    (params["alpha"] * par_seq + params["pmax"])

  r_squared <- 1 - sum(residuals(fit)^2) / sum((gpp_clean - mean(gpp_clean))^2)

  result <- list(
    alpha = params["alpha"],
    Pmax = params["pmax"],
    R2 = r_squared,
    par_data = par_clean,
    gpp_data = gpp_clean,
    par_pred = par_seq,
    gpp_pred = as.vector(gpp_pred),
    n_obs = length(par_clean)
  )

  if (!is.null(conf_int)) {
    result$alpha_se <- (conf_int["alpha", 2] - conf_int["alpha", 1]) / (2 * 1.96)
    result$Pmax_se <- (conf_int["pmax", 2] - conf_int["pmax", 1]) / (2 * 1.96)
  }

  return(result)
}


# Corrected wrapper function
fit_reco_gpp_from_nee <- function(temp, nee, par, par_threshold = 0.5) {
  # Complete workflow like MATLAB script

  cat("=== Step 1: Fitting Reco from nighttime NEE ===\n")
  reco_fit <- fit_reco_from_nee(temp, nee, par, par_threshold)

  cat("\n=== Step 2: Estimating GPP = Reco - NEE for ALL time points ===\n")
  estimates <- estimate_gpp_from_nee(temp, nee, reco_fit)

  cat("\n=== Step 3: Fitting GPP light response using ALL data ===\n")
  gpp_fit <- fit_gpp_lightresponse(par, estimates$gpp_est)  # No threshold!

  return(list(
    reco_fit = reco_fit,
    gpp_fit = gpp_fit,
    estimates = estimates
  ))
}
# Plotting functions (same as before)
plot_reco_response <- function(fit_result) {
  data_df <- data.frame(Temp = fit_result$temp_data, Reco = fit_result$reco_data)
  pred_df <- data.frame(Temp = fit_result$temp_pred, Reco = fit_result$reco_pred)

  r_ref_label <- sprintf("R[ref] == %.2f", fit_result$R_ref)
  e0_label <- sprintf("E[0] == %.0f", fit_result$E0)
  r2_label <- sprintf("R^2 == %.2f", fit_result$R2)

  temp_max <- max(data_df$Temp)
  reco_max <- max(data_df$Reco)

  ggplot() +
    geom_point(data = data_df, aes(x = Temp, y = Reco),
               alpha = 0.2, size = 1) +
    geom_line(data = pred_df, aes(x = Temp, y = Reco),
              color = "red", linewidth = 1) +
    annotate("text", x = temp_max * 0.95, y = reco_max * 0.95,
             label = r_ref_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = temp_max * 0.95, y = reco_max * 0.88,
             label = e0_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = temp_max * 0.95, y = reco_max * 0.81,
             label = r2_label, hjust = 1, size = 3.5, parse = TRUE) +
    labs(x = "Temperature (°C)",
         y = expression("Night-time NEE ("*mu*mol~CO[2]~m^-2~s^-1*")")) +
    theme_test()
}


# Improved plotting function with better spacing
plot_gpp_response <- function(fit_result) {
  data_df <- data.frame(PAR = fit_result$par_data, GPP = fit_result$gpp_data)
  pred_df <- data.frame(PAR = fit_result$par_pred, GPP = fit_result$gpp_pred)

  # Create labels with better formatting
  alpha_label <- sprintf("α == %.3f", fit_result$alpha)
  pmax_label <- sprintf("P[max] == %.2f", fit_result$Pmax)
  r2_label <- sprintf("R^2 == %.2f", fit_result$R2)

  # Calculate plot ranges with padding
  par_range <- range(data_df$PAR)
  gpp_range <- range(data_df$GPP)

  par_padding <- diff(par_range) * 0.05
  gpp_padding <- diff(gpp_range) * 0.05

  # Position text in top right with better spacing
  text_x <- par_range[2] - diff(par_range) * 0.05
  text_y_top <- gpp_range[2] - diff(gpp_range) * 0.05
  text_spacing <- diff(gpp_range) * 0.08

  ggplot() +
    geom_point(data = data_df, aes(x = PAR, y = GPP),
               alpha = 0.2, size = 0.8) +  # Smaller points, more transparent
    geom_line(data = pred_df, aes(x = PAR, y = GPP),
              color = "blue", linewidth = 1.2) +
    annotate("text", x = text_x, y = text_y_top,
             label = alpha_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = text_x, y = text_y_top - text_spacing,
             label = pmax_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = text_x, y = text_y_top - 2 * text_spacing,
             label = r2_label, hjust = 1, size = 3.5, parse = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = expression("PAR ("*mu*mol~m^-2~s^-1*")"),
         y = expression(GPP~"("*mu*mol~CO[2]~m^-2~s^-1*")")) +
    theme_test() +
    theme(
      plot.margin = margin(10, 15, 10, 10)
    )
}
