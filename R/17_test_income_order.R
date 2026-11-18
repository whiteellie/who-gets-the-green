source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

model_dir <- here::here("outputs", "models")
models <- readRDS(file.path(model_dir, "all_regions_fractional_logit_models.rds"))
reference <- read_csv(
  file.path(model_dir, "all_regions_fractional_logit_cv_predictions.csv"),
  show_col_types = FALSE
)
original_names <- names(models)[grepl("^M[1-5]:", names(models))]
stopifnot(length(original_names) == 5L)
final_reference <- filter(reference, model == original_names[[5]])
analysis <- model.frame(models[[original_names[[5]]]])
stopifnot(
  all(c("gpcd_log_z", "income_log_z", "precip_z", "eto_z", "developed_z", "space_z", "region") %in% names(analysis)),
  nrow(analysis) == nrow(final_reference),
  max(abs(100 * analysis$green_fraction - final_reference$observed_pct)) < 1e-8,
  !anyDuplicated(final_reference[c("pwsid", "region")]),
  all(final_reference$cv_fold %in% seq_len(5)),
  all(vapply(split(final_reference$cv_fold, final_reference$pwsid), function(folds) length(unique(folds)) == 1L, logical(1)))
)
analysis$cv_fold <- final_reference$cv_fold

reordered_formulas <- list(
  green_fraction ~ gpcd_log_z + region,
  green_fraction ~ gpcd_log_z + region + income_log_z,
  green_fraction ~ gpcd_log_z + region + income_log_z + precip_z + eto_z,
  green_fraction ~ gpcd_log_z + region + income_log_z + precip_z + eto_z + developed_z,
  formula(models[[original_names[[5]]]])
)
reordered_additions <- c("Water use + region", "Income", "Climate", "Development", "Potential landscape space")
original_additions <- c("Water use + region", "Climate", "Development", "Potential landscape space", "Income")

reordered_predictions <- bind_rows(lapply(seq_along(reordered_formulas), function(step) {
  predicted_pct <- rep(NA_real_, nrow(analysis))
  for (fold in seq_len(5)) {
    test_rows <- which(analysis$cv_fold == fold)
    fitted <- glm(
      reordered_formulas[[step]], family = quasibinomial(link = "logit"),
      data = analysis[analysis$cv_fold != fold, ]
    )
    stopifnot(fitted$converged, all(is.finite(coef(fitted))))
    predicted_pct[test_rows] <- 100 * predict(fitted, newdata = analysis[test_rows, ], type = "response")
  }
  transmute(
    final_reference, sequence = "Income second", step = .env$step,
    added_predictor = reordered_additions[[.env$step]],
    pwsid, region, cv_fold, observed_pct, predicted_pct = .env$predicted_pct
  )
}))
original_predictions <- reference |>
  filter(model %in% original_names) |>
  mutate(step = match(model, original_names)) |>
  transmute(
    sequence = "Income last", step, added_predictor = original_additions[step],
    pwsid, region, cv_fold, observed_pct, predicted_pct
  )
endpoint_check <- inner_join(
  filter(reordered_predictions, step %in% c(1L, 5L)),
  filter(original_predictions, step %in% c(1L, 5L)),
  by = c("step", "pwsid", "region", "cv_fold"), suffix = c("_reordered", "_original")
)
stopifnot(
  nrow(endpoint_check) == 2L * nrow(analysis),
  all(is.finite(reordered_predictions$predicted_pct)),
  max(abs(endpoint_check$predicted_pct_reordered - endpoint_check$predicted_pct_original)) < 1e-7
)

predictions <- bind_rows(original_predictions, reordered_predictions)
metrics <- predictions |>
  group_by(sequence, step, added_predictor) |>
  summarise(
    observations = n(), unique_suppliers = n_distinct(pwsid),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    .groups = "drop"
  ) |>
  arrange(sequence, step) |>
  group_by(sequence) |>
  mutate(
    incremental_r_squared = r_squared - lag(r_squared),
    rmse_reduction_percentage_points = lag(rmse_percentage_points) - rmse_percentage_points
  ) |>
  ungroup()

output_dir <- file.path(model_dir, "income_order_experiment")
figure_dir <- file.path(dir_figures, "archived", "income_order_experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(metrics, file.path(output_dir, "cv_metrics.csv"))
write_csv(predictions, file.path(output_dir, "cv_predictions.csv"))
write_csv(filter(metrics, added_predictor == "Income"), file.path(output_dir, "income_increment.csv"))

comparison_plot <- ggplot(metrics, aes(step, r_squared, color = sequence, group = sequence)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sprintf("%.3f", r_squared)), vjust = -0.8, show.legend = FALSE, check_overlap = TRUE) +
  facet_wrap(~sequence) +
  scale_x_continuous(breaks = seq_len(5)) +
  scale_color_manual(values = c("Income last" = "#176b3a", "Income second" = "#b85c26")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Does income contribute more when added earlier?",
    subtitle = paste(nrow(analysis), "records; identical sample and saved PWS-grouped five-fold splits"),
    x = "Cumulative model step", y = "Out-of-fold R-squared",
    caption = paste(
      "Income last: water use + region -> climate -> development -> landscape space -> income.",
      "Income second: water use + region -> income -> climate -> development -> landscape space.",
      "Sequential gains depend on order; the final predictor set is unchanged.", sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", plot.caption = element_text(hjust = 0))
ggsave(file.path(figure_dir, "income_order_cv_comparison.png"), comparison_plot, width = 10, height = 6, dpi = 300)
print(metrics, n = Inf, width = Inf)
cat("Baseline and final predictions reproduce retained results within 1e-7 percentage points.\n")