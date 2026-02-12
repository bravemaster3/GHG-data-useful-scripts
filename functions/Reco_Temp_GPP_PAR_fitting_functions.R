library(ggplot2)

# Function 1: Lloyd & Taylor (1994) temperature response for Reco https://doi.org/10.2307/2389824
# Reco = R10 * exp(E0 * (1/(Tref - T0) - 1/(Tsoil - T0)))
fit_reco_lloyd_taylor <- function(temp, reco) {
  # temp: temperature (°C) - can be air temp, soil temp, etc.
  # reco: ecosystem respiration (µmol CO2 m-2 s-1)

  # Remove missing values
  valid <- !is.na(temp) & !is.na(reco) #& reco > 0
  temp_clean <- temp[valid]
  reco_clean <- reco[valid]

  # Convert to Kelvin
  temp_K <- temp_clean + 273.15
  tref <- 283.15  # 10°C in Kelvin
  t0 <- 227.13    # -46.02°C in Kelvin (Lloyd & Taylor 1994)

  # Fit model
  fit <- nls(reco_clean ~ r10 * exp(e0 * (1/(tref - t0) - 1/(temp_K - t0))),
             start = list(r10 = median(reco_clean), e0 = 400),
             control = nls.control(maxiter = 200))

  # Get parameters and confidence intervals
  params <- coef(fit)
  conf_int <- confint(fit)

  # Predictions for smooth line
  temp_seq <- seq(min(temp_clean), max(temp_clean), length.out = 100)
  temp_seq_K <- temp_seq + 273.15
  reco_pred <- params["r10"] * exp(params["e0"] *
                                     (1/(tref - t0) - 1/(temp_seq_K - t0)))

  return(list(
    R10 = params["r10"],
    R10_se = (conf_int["r10", 2] - conf_int["r10", 1]) / (2 * 1.96),
    E0 = params["e0"],
    E0_se = (conf_int["e0", 2] - conf_int["e0", 1]) / (2 * 1.96),
    R2 = 1 - sum(residuals(fit)^2) / sum((reco_clean - mean(reco_clean))^2),
    temp_data = temp_clean,
    reco_data = reco_clean,
    temp_pred = temp_seq,
    reco_pred = as.vector(reco_pred),
    n_obs = length(temp_clean)
  ))
}

# Function 3: Hyperbolic light response for GPP [Michaelis-Menten rectangular hyperbola]
# GPP = -α * PAR * Pmax / (α * PAR + Pmax)

fit_gpp_hyperbolic <- function(par, gpp, par_min_threshold=0) {
  # par: photosynthetically active radiation (µmol m-2 s-1)
  # gpp: gross primary productivity (µmol CO2 m-2 s-1, NEGATIVE values)

  # Remove missing values and nighttime
  valid <- !is.na(par) & !is.na(gpp) & par > par_min_threshold & gpp < 0
  par_clean <- par[valid]
  gpp_clean <- abs(gpp[valid])  # Work with positive values

  # Simple starting values
  alpha_init <- 0.05
  pmax_init <- max(gpp_clean, na.rm = TRUE)

  # Fit model with EXPLICIT FORMULA
  # GPP = (alpha * PAR * Pmax) / (alpha * PAR + Pmax)
  fit <- nls(gpp_clean ~ (alpha * par_clean * pmax) / (alpha * par_clean + pmax),
             start = list(alpha = alpha_init, pmax = pmax_init),
             control = nls.control(maxiter = 1000, warnOnly = TRUE))

  # Get parameters
  params <- coef(fit)
  conf_int <- tryCatch(confint(fit), error = function(e) NULL)

  # Predictions
  par_seq <- seq(0, max(par_clean), length.out = 100)
  gpp_pred <- (params["alpha"] * par_seq * params["pmax"]) /
    (params["alpha"] * par_seq + params["pmax"])

  result <- list(
    alpha = params["alpha"],
    Pmax = params["pmax"],
    R2 = 1 - sum(residuals(fit)^2) / sum((gpp_clean - mean(gpp_clean))^2),
    par_data = par_clean,
    gpp_data = -gpp_clean,  # Return as negative
    par_pred = par_seq,
    gpp_pred = -gpp_pred,   # Return as negative
    n_obs = length(par_clean)
  )

  if (!is.null(conf_int)) {
    result$alpha_se <- (conf_int["alpha", 2] - conf_int["alpha", 1]) / (2 * 1.96)
    result$Pmax_se <- (conf_int["pmax", 2] - conf_int["pmax", 1]) / (2 * 1.96)
  }

  return(result)
}


# Plotting functions
plot_reco_response <- function(fit_result) {
  data_df <- data.frame(
    Temp = fit_result$temp_data,
    Reco = fit_result$reco_data
  )

  pred_df <- data.frame(
    Temp = fit_result$temp_pred,
    Reco = fit_result$reco_pred
  )

  # Create labels
  if (!is.null(fit_result$R10_se)) {
    r10_label <- sprintf("R[10] == %.2f", fit_result$R10)
  } else {
    r10_label <- sprintf("R[10] == %.2f", fit_result$R10)
  }

  if (!is.null(fit_result$E0_se)) {
    e0_label <- sprintf("E[0] == %.0f", fit_result$E0)
  } else {
    e0_label <- sprintf("E[0] == %.0f", fit_result$E0)
  }

  r2_label <- sprintf("R^2 == %.2f", fit_result$R2)

  # Get plot ranges for text positioning
  temp_max <- max(data_df$Temp)
  reco_max <- max(data_df$Reco)

  ggplot() +
    geom_point(data = data_df, aes(x = Temp, y = Reco),
               alpha = 0.2, size = 1) +
    geom_line(data = pred_df, aes(x = Temp, y = Reco),
              color = "red", linewidth = 1) +
    annotate("text", x = temp_max * 0.15, y = reco_max * 0.95,
             label = r10_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = temp_max * 0.15, y = reco_max * 0.88,
             label = e0_label, hjust = 1, size = 3.5, parse = TRUE) +
    annotate("text", x = temp_max * 0.15, y = reco_max * 0.81,
             label = r2_label, hjust = 1, size = 3.5, parse = TRUE) +
    labs(x = "Temperature (°C)",
         y = expression(R[eco]~"("*mu*mol~CO[2]~m^-2~s^-1*")")) +
    theme_test()
}


plot_gpp_response <- function(fit_result) {
  # Create data frames for plotting
  data_df <- data.frame(
    PAR = fit_result$par_data,
    GPP = fit_result$gpp_data
  )

  pred_df <- data.frame(
    PAR = fit_result$par_pred,
    GPP = fit_result$gpp_pred
  )

  # Create label text
  if (!is.null(fit_result$alpha_se)) {
    alpha_label <- sprintf("α = %.3f ± %.3f", fit_result$alpha, fit_result$alpha_se)
  } else {
    alpha_label <- sprintf("α = %.3f", fit_result$alpha)
  }

  if (!is.null(fit_result$Pmax_se)) {
    pmax_label <- sprintf("P[max] = %.2f ± %.2f", fit_result$Pmax, fit_result$Pmax_se)
  } else {
    pmax_label <- sprintf("P[max] = %.2f", fit_result$Pmax)
  }

  r2_label <- sprintf("R² = %.2f", fit_result$R2)

  # Combine labels
  labels_text <- paste(alpha_label, pmax_label, r2_label, sep = "\n")

  # Create plot
  ggplot() +
    geom_point(data = data_df, aes(x = PAR, y = GPP),
               alpha = 0.2, size = 1) +
    geom_line(data = pred_df, aes(x = PAR, y = GPP),
              color = "blue", linewidth = 1) +
    annotate("text", x = 1900, y = -0.5, label = labels_text,
             hjust = 1, vjust = 1, size = 3.5) +
    scale_x_continuous(limits = c(0, 2000)) +
    scale_y_continuous(limits = c(-8, 0)) +
    labs(x = expression("PAR ("*mu*mol~m^-2~s^-1*")"),
         y = expression(GPP~"("*mu*mol~CO[2]~m^-2~s^-1*")")) +
    theme_test() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 11))

}


# Example usage:
# reco_fit <- fit_reco_lloyd_taylor(tsoil = df$Ts10, reco = df$Reco)
# gpp_fit <- fit_gpp_hyperbolic(par = df$PAR, gpp = df$GPP)
#
# par(mfrow = c(1, 2))
# plot_reco_response(reco_fit)
# plot_gpp_response(gpp_fit)
