library(tidyverse)
library(dplyr)
library(ggplot2)
library(zoo)

# ── Settings ──────────────────────────────────────────────────────────────────

roll_days <- 7  # change rolling window size here (number of days averaged)

site_colors <- c(
  "SE-Deg" = "#E69F00",
  "SE-HfM" = "#56B4E9",
  "SE-Hmr" = "#009E73",
  "SE-Srj" = "#CC79A7"
)

# ── Custom x-axis formatter: J F M A M J J A S O N D ────────────────────────

dte_formatter <- function(x) {
  mth <- toupper(substr(format(x, "%b"), 1, 1))
  ifelse(format(x, "%m") == "01", paste0(mth, "\n", format(x, "%Y")), mth)
}

# ── Load & compute rolling means ─────────────────────────────────────────────

daily <- read_csv("data/daily_data/all_sites_daily.csv",
                  col_types = cols(date = col_date(), site = col_character(),
                                   .default = col_double())) %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    tair_roll       = rollmean(tair,       k = roll_days, fill = NA, align = "center"),
    tsoil_10cm_roll = rollmean(tsoil_10cm, k = roll_days, fill = NA, align = "center"),
    WTD_roll        = rollmean(WTD,        k = roll_days, fill = NA, align = "center"),
    PARin_roll      = rollmean(PARin,      k = roll_days, fill = NA, align = "center")
  ) %>%
  ungroup()

# ── Pivot to long ─────────────────────────────────────────────────────────────

var_labels <- c(
  tair       = '"Air temp. (°C)"',
  tsoil_10cm = '"Soil temp. 10 cm (°C)"',
  WTD        = '"WTD (cm)"',
  PARin      = 'PAR["in"]~(μmol~m^{-2}~s^{-1})'
)

long <- daily %>%
  pivot_longer(c(tair, tsoil_10cm, WTD, PARin),
               names_to = "variable", values_to = "daily_mean") %>%
  pivot_longer(c(tair_roll, tsoil_10cm_roll, WTD_roll, PARin_roll),
               names_to = "variable_roll", values_to = "roll") %>%
  filter(variable_roll == paste0(variable, "_roll")) %>%
  mutate(variable = factor(variable, levels = names(var_labels)))

# ── Plot ──────────────────────────────────────────────────────────────────────

ggplot(long, aes(x = date, color = site)) +
  geom_point(aes(y = daily_mean), alpha = 0.5, size = 0.6, shape = 16) +
  geom_line(aes(y = roll), linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y",
             labeller = as_labeller(var_labels, label_parsed)) +
  scale_color_manual(values = site_colors) +
  scale_x_date(date_breaks = "1 month", labels = dte_formatter, expand = c(0, 0)) +
  labs(x = NULL, y = "Daily mean", color = "Site") +
  theme_test() +
  theme(
    legend.position = "top",
    strip.text      = element_text(size = 9, hjust = 0.5),
    axis.text.x     = element_text(size = 8)
  )

ggsave("figures/biomet_daily_roll.png", width = 10, height = 10, dpi = 300)
message("Plot saved to figures/biomet_daily_roll.png")
