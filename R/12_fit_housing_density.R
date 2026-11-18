source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("sandwich", quietly = TRUE)) stop("Install the 'sandwich' package first.")
fold_path <- here::here("outputs", "models", "all_regions_fractional_logit_cv_predictions.csv")
if (!file.exists(fold_path)) stop("Run the retained R/09_fit_models.R analysis first.")
retained_cv <- read_csv(fold_path, show_col_types = FALSE)
retained_primary <- retained_cv |>
  filter(grepl("^M1:", model)) |>
  select(pwsid, region, cv_fold)
stopifnot(nrow(retained_primary) > 0L, !anyDuplicated(retained_primary[c("pwsid", "region")]))
stopifnot(all(vapply(split(retained_primary$cv_fold, retained_primary$pwsid), function(folds) n_distinct(folds) == 1L, logical(1))))

acs <- sf::read_sf(file.path(dir_raw, "acs_bg.gpkg")) |>
  sf::st_drop_geometry() |>
  select(GEOID, housing_units)
acs$GEOID <- as.character(acs$GEOID)
stopifnot(!anyDuplicated(acs$GEOID))
region_dirs <- sort(list.dirs(dir_processed, recursive = FALSE, full.names = TRUE))
region_dirs <- region_dirs[grepl("^hr_", basename(region_dirs))]
stopifnot(length(region_dirs) == 10L)

density_components <- bind_rows(lapply(region_dirs, function(region_dir) {
  scope <- basename(region_dir)
  overlaps <- read_csv(
    file.path(region_dir, "block_group_pws.csv"),
    col_types = cols(GEOID = col_character(), .default = col_guess()),
    col_select = all_of(c("GEOID", "pwsid", "overlap_area_m2", "bg_area", "w"))
  ) |>
    left_join(acs, by = "GEOID")
  stopifnot(all(is.finite(overlaps$overlap_area_m2)), all(overlaps$overlap_area_m2 > 0))
  stopifnot(all(is.finite(overlaps$bg_area)), all(overlaps$bg_area > 0), all(overlaps$w > 0))
  stopifnot(all(overlaps$overlap_area_m2 / overlaps$bg_area <= 1.0001))
  overlaps |>
    mutate(
      region = tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope))),
      housing_valid = is.finite(housing_units) & housing_units >= 0,
      allocated_housing_units = if_else(housing_valid, housing_units * overlap_area_m2 / bg_area, NA_real_)
    ) |>
    select(region, pwsid, GEOID, housing_units, housing_valid, overlap_area_m2, bg_area, allocated_housing_units)
}))
density <- density_components |>
  group_by(region, pwsid) |>
  summarise(
    sampled_footprint_area_km2 = sum(overlap_area_m2) / 1e6,
    housing_coverage_pct = 100 * sum(overlap_area_m2[housing_valid]) / sum(overlap_area_m2),
    housing_complete = all(housing_valid),
    estimated_housing_units = if (all(housing_valid)) sum(allocated_housing_units) else NA_real_,
    .groups = "drop"
  ) |>
  mutate(housing_units_per_km2 = estimated_housing_units / sampled_footprint_area_km2)

supplier_rows <- bind_rows(lapply(region_dirs, function(region_dir) {
  scope <- basename(region_dir)
  read_csv(file.path(region_dir, "pws_green_water_income.csv"), show_col_types = FALSE) |>
    mutate(region = tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope))))
}))
stopifnot(!anyDuplicated(supplier_rows[c("pwsid", "region")]))
candidates <- retained_primary |>
  left_join(supplier_rows, by = c("pwsid", "region")) |>
  left_join(density, by = c("pwsid", "region")) |>
  mutate(density_eligible = housing_complete %in% TRUE & is.finite(housing_units_per_km2) & housing_units_per_km2 >= 0)
analysis <- candidates |>
  filter(density_eligible) |>
  mutate(
    green_fraction = pct_green_sampled / 100,
    log_outdoor = log1p(outdoor_gpcd),
    log_income = log(income_proxy),
    log_density = log1p(housing_units_per_km2),
    region = factor(region)
  )
continuous <- c("log_outdoor", "precip_gs_mm", "eto_gs_mm", "developed_fraction_sampled", "log_income", "log_density")
stopifnot(nrow(analysis) > 0L, all(between(analysis$outdoor_gpcd, 0, 132)))
stopifnot(all(between(analysis$green_fraction, 0, 1)))
stopifnot(all(vapply(analysis[continuous], function(values) all(is.finite(values)) && sd(values) > 0, logical(1))))

formulas <- list(
  "M3: development, no income" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + region,
  "M3 + housing density" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_density + region,
  "M5: income, no density" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + region,
  "M5 + housing density" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + log_density + region
)
fit_model <- function(data, formula) {
  fit <- glm(formula, family = quasibinomial(link = "logit"), data = data)
  if (!fit$converged || any(!is.finite(coef(fit)))) stop("Model did not converge with a full-rank design.")
  fit
}
models <- lapply(formulas, function(formula) fit_model(analysis, formula))
cv_predictions <- bind_rows(lapply(names(formulas), function(model_name) {
  predicted <- rep(NA_real_, nrow(analysis))
  for (fold in seq_len(5)) {
    test <- which(analysis$cv_fold == fold)
    train <- analysis[analysis$cv_fold != fold, ]
    stopifnot(n_distinct(train$region) == n_distinct(analysis$region))
    fit <- fit_model(train, formulas[[model_name]])
    predicted[test] <- 100 * predict(fit, newdata = analysis[test, ], type = "response")
  }
  data.frame(
    model = model_name, pwsid = analysis$pwsid, region = analysis$region,
    fold = analysis$cv_fold, observed_pct = analysis$pct_green_sampled, predicted_pct = predicted
  )
}))
stopifnot(all(is.finite(cv_predictions$predicted_pct)))
cv_predictions$model <- factor(cv_predictions$model, levels = names(formulas))
cv_metrics <- cv_predictions |>
  group_by(model) |>
  summarise(
    observations = n(), pws_ids = n_distinct(pwsid),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    mae_percentage_points = mean(abs(observed_pct - predicted_pct)), .groups = "drop"
  )
comparison <- data.frame(
  reference_model = c("M3: development, no income", "M5: income, no density"),
  density_model = c("M3 + housing density", "M5 + housing density")
)
comparison$reference_r_squared <- cv_metrics$r_squared[match(comparison$reference_model, cv_metrics$model)]
comparison$density_r_squared <- cv_metrics$r_squared[match(comparison$density_model, cv_metrics$model)]
comparison$delta_r_squared <- comparison$density_r_squared - comparison$reference_r_squared
comparison$reference_rmse <- cv_metrics$rmse_percentage_points[match(comparison$reference_model, cv_metrics$model)]
comparison$density_rmse <- cv_metrics$rmse_percentage_points[match(comparison$density_model, cv_metrics$model)]

log_density_sd <- sd(analysis$log_density)
density_coefficients <- bind_rows(lapply(comparison$density_model, function(model_name) {
  fit <- models[[model_name]]
  covariance <- sandwich::vcovCL(fit, cluster = analysis$pwsid, type = "HC1")
  estimate <- coef(fit)[["log_density"]] * log_density_sd
  se <- sqrt(covariance["log_density", "log_density"]) * log_density_sd
  data.frame(
    model = model_name, log_density_sd = log_density_sd,
    estimate_log_odds_per_sd = estimate, cluster_robust_se = se,
    odds_ratio_per_sd = exp(estimate), ci_low = exp(estimate - 1.96 * se),
    ci_high = exp(estimate + 1.96 * se), p_value = 2 * pnorm(-abs(estimate / se))
  )
}))

full_model <- models[["M5 + housing density"]]
beta <- coef(full_model)
full_vcov <- sandwich::vcovCL(full_model, cluster = analysis$pwsid, type = "HC1")
density_grid <- expm1(seq(
  quantile(analysis$log_density, 0.05), quantile(analysis$log_density, 0.95), length.out = 80
))
adjusted <- bind_rows(lapply(density_grid, function(value) {
  new_data <- analysis
  new_data$log_density <- log1p(value)
  design <- model.matrix(delete.response(terms(full_model)), new_data)
  design <- design[, names(beta), drop = FALSE]
  response <- plogis(drop(design %*% beta))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  se <- sqrt(drop(t(gradient) %*% full_vcov %*% gradient))
  data.frame(
    housing_units_per_km2 = value, fitted_pct_green = 100 * mean(response),
    ci_low = 100 * max(0, mean(response) - 1.96 * se),
    ci_high = 100 * min(1, mean(response) + 1.96 * se)
  )
}))

output_dir <- here::here("outputs", "models", "housing_density_experiment")
figure_dir <- file.path(dir_figures, "archived", "housing_density_experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(density_components, file.path(output_dir, "density_overlap_components.csv"))
write_csv(density, file.path(output_dir, "supplier_region_density.csv"))
write_csv(candidates, file.path(output_dir, "sample_eligibility.csv"))
write_csv(analysis, file.path(output_dir, "model_sample.csv"))
write_csv(cv_predictions, file.path(output_dir, "cv_predictions.csv"))
write_csv(cv_metrics, file.path(output_dir, "cv_metrics.csv"))
write_csv(comparison, file.path(output_dir, "incremental_cv_comparison.csv"))
write_csv(density_coefficients, file.path(output_dir, "density_coefficients.csv"))
write_csv(adjusted, file.path(output_dir, "adjusted_density_predictions.csv"))
write_csv(count(analysis, region, name = "observations"), file.path(output_dir, "sample_by_region.csv"))
write_csv(
  data.frame(
    original_primary_records = nrow(retained_primary), density_complete_records = nrow(analysis),
    excluded_for_density = nrow(candidates) - nrow(analysis),
    source = "ACS B25001 housing units already in acs_bg.gpkg",
    aggregation = "sum(housing_units * overlap_area_m2 / original_bg_area_m2) / sum(overlap_area_m2) * 1e6",
    footprint = "Existing sampled positive-population supplier-region overlaps; gross polygon area, including water/nonresidential land",
    assumption = "Housing units uniformly distributed within each original block group; no missing housing counts filled with zero",
    predictor = "log1p housing units per square kilometer; coefficients reported per SD of transformed density",
    validation = "Exact retained-model PWS fold assignments; all four models use the same density-complete rows"
  ), file.path(output_dir, "analysis_definition.csv")
)
saveRDS(models, file.path(output_dir, "models.rds"))

metric_labels <- cv_metrics |>
  mutate(label = paste0("R2 = ", sprintf("%.3f", r_squared), "\nRMSE = ", sprintf("%.1f", rmse_percentage_points), " pp"))
prediction_plot <- ggplot(cv_predictions, aes(observed_pct, predicted_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_point(color = "#2c7fb8", alpha = 0.25, size = 1.2) +
  geom_text(data = metric_labels, aes(x = 3, y = 97, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  coord_equal() +
  labs(
    x = "Observed developed-land greenness (%)", y = "Out-of-fold predicted greenness (%)",
    title = "Does housing density improve prediction?",
    subtitle = paste(nrow(analysis), "matched records; <=132 GPCD; retained PWS-grouped folds"),
    caption = "Density uses allocated housing units divided by total sampled area, not an unweighted average of block-group densities.\nAll panels use identical records and folds; the retained primary model is unchanged."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "predicted_vs_observed.png"), prediction_plot, width = 9, height = 9, dpi = 180)

effect_plot <- ggplot() +
  geom_point(data = analysis, aes(log1p(housing_units_per_km2), pct_green_sampled), alpha = 0.15, color = "#777777", size = 1) +
  geom_ribbon(data = adjusted, aes(log1p(housing_units_per_km2), ymin = ci_low, ymax = ci_high), fill = "#2c7fb8", alpha = 0.18) +
  geom_line(data = adjusted, aes(log1p(housing_units_per_km2), fitted_pct_green), color = "#2c7fb8", linewidth = 1) +
  scale_x_continuous(breaks = log1p(c(0, 1, 10, 100, 1000, 10000)), labels = scales::label_comma()(c(0, 1, 10, 100, 1000, 10000))) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Housing units per square kilometer of sampled footprint (log1p-scaled axis)",
    y = "Greenness within developed land (%)",
    title = "Adjusted greenness across housing density",
    subtitle = "M5 + density: marginal predictions averaging over other observed covariates and regions",
    caption = "Curve shown over the 5th-95th percentiles; shading is a 95% PWS-clustered interval.\nGross density is development context, not irrigable area; associations are not causal."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "adjusted_greenness_vs_housing_density.png"), effect_plot, width = 10, height = 7, dpi = 180)

cat("Original retained records:", nrow(retained_primary), "density-complete records:", nrow(analysis), "\n")
print(cv_metrics, width = Inf)
print(comparison, row.names = FALSE)
print(density_coefficients, row.names = FALSE)
cat("Saved housing-density experiment to", output_dir, "and", figure_dir, "\n")