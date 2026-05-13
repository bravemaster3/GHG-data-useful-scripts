library(data.table)
library(dplyr)
library(ggplot2)

abisko_df <- fread("data/Karins/abisko_df.csv")
abisko_raw <- fread("data/Karins/abisko_raw.csv")

fit_reco_lloyd_taylor <- function(temp, reco) {

  temp <- unname(as.numeric(temp))
  reco <- unname(as.numeric(reco))

  valid <- !is.na(temp) & !is.na(reco) & reco > 0
  temp_clean <- temp[valid]
  reco_clean <- reco[valid]

  temp_K <- temp_clean + 273.15
  tref <- 283.15
  t0 <- 227.13

  fit <- nls(
    reco_clean ~ r10 * exp(e0 * (1/(tref - t0) - 1/(temp_K - t0))),
    start = list(r10 = median(reco_clean, na.rm = TRUE), e0 = 400),
    control = nls.control(maxiter = 1000, warnOnly = TRUE)
  )

  params <- coef(fit)

  temp_seq <- seq(min(temp_clean), max(temp_clean), length.out = 100)
  temp_seq_K <- temp_seq + 273.15
  reco_pred <- params["r10"] * exp(params["e0"] *
                                     (1/(tref - t0) - 1/(temp_seq_K - t0)))

  result <- list(
    R10 = params["r10"],
    E0 = params["e0"],
    R2 = 1 - sum(residuals(fit)^2) / sum((reco_clean - mean(reco_clean))^2),
    temp_data = temp_clean,
    reco_data = reco_clean,
    temp_pred = temp_seq,
    reco_pred = as.vector(reco_pred),
    n_obs = length(temp_clean),
    model = fit
  )

  return(result)
}

# abisko_df <- all_sites %>%
#   filter(site == "Abisko")

abisko_reco_df <- abisko_df %>%
  filter(
    !is.na(temp),
    !is.na(reco),
    par < 10,
    temp > - 10 #### filtered low temperatures (< −10 °C)
  )

abisko_reco_fit <- fit_reco_lloyd_taylor(
  temp = abisko_reco_df$temp,
  reco = abisko_reco_df$nee
)

plot(nee~reco, data=abisko_reco_df)
abline(0,1, col=2)

plot_reco_response <- function(fit_result) {

  data_df <- data.frame(
    Temp = fit_result$temp_data,
    Reco = fit_result$reco_data
  )

  pred_df <- data.frame(
    Temp = fit_result$temp_pred,
    Reco = fit_result$reco_pred
  )

  ggplot() +
    geom_point(data = data_df, aes(x = Temp, y = Reco),
               alpha = 0.2, size = 1) +
    geom_line(data = pred_df, aes(x = Temp, y = Reco),
              color = "red", linewidth = 1) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(
      x = "Temperature (°C)",
      y = expression(R[eco]~"("*mu*mol~CO[2]~m^-2~s^-1*")")
    ) +
    theme_test()
}

g <- plot_reco_response(abisko_reco_fit)
g
# g + ylim(0,5)
