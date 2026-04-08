library(data.table)
library(ggplot2)
library(dplyr)

height <- "ws100"
input_dir <- file.path("results/univariate_garch", height)

results_dt <- readRDS(file.path(input_dir, "garch_results.rds"))

# -----------------------------
# GARCH plot
# -----------------------------
garch_results <- results_dt[model == "sGARCH"]

garch_results[, sum_ab := alpha1 + beta1]

gr <- garch_results %>%
  filter(!is.na(alpha1), !is.na(beta1)) %>%
  mutate(
    persist = ifelse(sum_ab < 1, "<1", "≥1")
  )

p_garch <- ggplot(gr, aes(alpha1, beta1, colour = persist)) +
  geom_point(size = 2.5) +
  theme_minimal() +
  labs(title = "GARCH coefficients", x = "alpha", y = "beta")

ggsave("garch_plot.pdf", p_garch, width = 6, height = 4)

# -----------------------------
# EGARCH plot
# -----------------------------
eg <- results_dt[model == "eGARCH"]

eg_df <- eg %>%
  filter(!is.na(alpha1), !is.na(beta1))

p_egarch <- ggplot(eg_df, aes(alpha1, beta1)) +
  geom_point() +
  theme_minimal() +
  labs(title = "EGARCH coefficients")

ggsave("egarch_plot.pdf", p_egarch, width = 6, height = 4)