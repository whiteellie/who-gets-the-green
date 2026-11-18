source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("sandwich", quietly = TRUE)) stop("Install the 'sandwich' package first.")
cv_path <- here::here("outputs", "models", "all_regions_fractional_logit_cv_predictions.csv")
if (!file.exists(cv_path)) stop("Run retained R/09_fit_models.R first.")
retained <- read_csv(cv_path, show_col_types = FALSE) |>
  filter(grepl("^M1:", model)) |>
  select(pwsid, region, cv_fold)
stopifnot(!anyDuplicated(retained[c("pwsid", "region")]))
acs <- sf::read_sf(file.path(dir_raw, "acs_bg.gpkg")) |>
  sf::st_drop_geometry() |>
  select(GEOID, housing_units)
acs$GEOID <- as.character(acs$GEOID)
stopifnot(!anyDuplicated(acs$GEOID))
regions <- sort(list.dirs(dir_processed, recursive = FALSE, full.names = TRUE))
regions <- regions[grepl("^hr_", basename(regions))]
stopifnot(length(regions) == 10L)

components <- bind_rows(lapply(regions, function(region_dir) {
  scope <- basename(region_dir)
  region_name <- tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope)))
  suppliers <- filter(retained, region == region_name)$pwsid
  checkpoint <- file.path(dir_raw, paste0("landscape_space_bg_", scope, ".csv"))
  temporary <- sub("[.]csv$", ".tmp", checkpoint)
  if (!file.exists(checkpoint) && !file.exists(temporary)) stop("Extract NLCD landscape space for ", region_name, " first.")
  nlcd <- if (file.exists(checkpoint)) read_csv(checkpoint, col_types = cols(GEOID = "c", .default = col_guess())) else NULL
  if (file.exists(temporary)) {
    recovered <- read_csv(temporary, col_types = cols(GEOID = "c", .default = col_guess()))
    if (!is.null(nlcd) && !all(nlcd$GEOID %in% recovered$GEOID)) stop("Incomplete temporary checkpoint: ", temporary)
    nlcd <- recovered
  }
  stopifnot(!anyDuplicated(nlcd$GEOID))
  overlap <- read_csv(
    file.path(region_dir, "block_group_pws.csv"),
    col_types = cols(GEOID = "c", .default = col_guess()),
    col_select = all_of(c("GEOID", "pwsid", "pct_green_of_developed", "overlap_area_m2", "bg_area", "w"))
  ) |>
    filter(pwsid %in% suppliers) |>
    left_join(nlcd, by = "GEOID") |>
    left_join(acs, by = "GEOID") |>
    mutate(
      region = region_name,
      green_valid = is.finite(pct_green_of_developed),
      space_valid = status %in% "ok" & is.finite(potential_space_fraction) &
        between(potential_space_fraction, 0, 1) & impervious_coverage >= 0.99,
      housing_valid = is.finite(housing_units) & housing_units >= 0,
      allocated_housing_units = if_else(housing_valid, housing_units * overlap_area_m2 / bg_area, NA_real_)
    )
  stopifnot(all(is.finite(overlap$w)), all(overlap$w > 0), all(overlap$overlap_area_m2 > 0), all(overlap$bg_area > 0))
  overlap
}))
proxy <- components |>
  group_by(pwsid, region) |>
  summarise(
    green_contributing_population_weight = sum(w[green_valid]),
    proxy_population_coverage = if (sum(w[green_valid]) > 0) sum(w[green_valid & space_valid]) / sum(w[green_valid]) else NA_real_,
    potential_space_fraction_sampled = if (any(green_valid & space_valid))
      weighted.mean(potential_space_fraction[green_valid & space_valid], w[green_valid & space_valid]) else NA_real_,
    sampled_footprint_area_km2 = sum(overlap_area_m2) / 1e6,
    housing_complete = all(housing_valid),
    housing_units_per_km2 = if (all(housing_valid)) sum(allocated_housing_units) / (sum(overlap_area_m2) / 1e6) else NA_real_,
    .groups = "drop"
  )
supplier_rows <- bind_rows(lapply(regions, function(region_dir) {
  scope <- basename(region_dir)
  read_csv(file.path(region_dir, "pws_green_water_income.csv"), show_col_types = FALSE) |>
    mutate(region = tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope))))
}))
candidates <- retained |>
  left_join(supplier_rows, by = c("pwsid", "region")) |>
  left_join(proxy, by = c("pwsid", "region")) |>
  mutate(
    proxy_eligible = is.finite(proxy_population_coverage) & proxy_population_coverage >= 0.99 &
      is.finite(potential_space_fraction_sampled) & housing_complete %in% TRUE &
      is.finite(housing_units_per_km2) & housing_units_per_km2 >= 0
  )
analysis <- candidates |>
  filter(proxy_eligible) |>
  mutate(
    green_fraction = pct_green_sampled / 100,
    log_outdoor = log1p(outdoor_gpcd),
    log_income = log(income_proxy),
    log_density = log1p(housing_units_per_km2),
    region = factor(region)
  )
predictors <- c("log_outdoor", "precip_gs_mm", "eto_gs_mm", "developed_fraction_sampled", "log_income", "log_density", "potential_space_fraction_sampled")
stopifnot(nrow(analysis) > 0L, all(between(analysis$outdoor_gpcd, 0, 132)))
stopifnot(all(between(analysis$green_fraction, 0, 1)))
stopifnot(all(vapply(analysis[predictors], function(values) all(is.finite(values)) && sd(values) > 0, logical(1))))

formulas <- list(
  "M3: development" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + region,
  "M3 + landscape space" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + potential_space_fraction_sampled + region,
  "M5: income" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + region,
  "M5 + landscape space" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + potential_space_fraction_sampled + region,
  "M5 + housing density" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + log_density + region,
  "M5 + space + density" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + log_density + potential_space_fraction_sampled + region
)
fit_model <- function(data, formula) {
  fit <- glm(formula, family = quasibinomial(link = "logit"), data = data)
  if (!fit$converged || any(!is.finite(coef(fit)))) stop("Model failed to converge with a full-rank design.")
  fit
}
models <- lapply(formulas, function(formula) fit_model(analysis, formula))
cv_predictions <- bind_rows(lapply(names(formulas), function(model_name) {
  predictions <- rep(NA_real_, nrow(analysis))
  for (fold in seq_len(5)) {
    test <- which(analysis$cv_fold == fold)
    train <- analysis[analysis$cv_fold != fold, ]
    stopifnot(n_distinct(train$region) == n_distinct(analysis$region))
    model <- fit_model(train, formulas[[model_name]])
    predictions[test] <- 100 * predict(model, newdata = analysis[test, ], type = "response")
  }
  data.frame(
    model = model_name, pwsid = analysis$pwsid, region = analysis$region,
    fold = analysis$cv_fold, observed_pct = analysis$pct_green_sampled, predicted_pct = predictions
  )
}))
stopifnot(all(is.finite(cv_predictions$predicted_pct)))
cv_predictions$model <- factor(cv_predictions$model, levels = names(formulas))
metrics <- cv_predictions |>
  group_by(model) |>
  summarise(
    observations = n(), pws_ids = n_distinct(pwsid),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    mae_percentage_points = mean(abs(observed_pct - predicted_pct)), .groups = "drop"
  )
comparisons <- data.frame(
  reference_model = c("M3: development", "M5: income", "M5 + housing density"),
  proxy_model = c("M3 + landscape space", "M5 + landscape space", "M5 + space + density")
)
comparisons$reference_r_squared <- metrics$r_squared[match(comparisons$reference_model, metrics$model)]
comparisons$proxy_r_squared <- metrics$r_squared[match(comparisons$proxy_model, metrics$model)]
comparisons$delta_r_squared <- comparisons$proxy_r_squared - comparisons$reference_r_squared
comparisons$reference_rmse <- metrics$rmse_percentage_points[match(comparisons$reference_model, metrics$model)]
comparisons$proxy_rmse <- metrics$rmse_percentage_points[match(comparisons$proxy_model, metrics$model)]

space_sd <- sd(analysis$potential_space_fraction_sampled)
coefficients <- bind_rows(lapply(comparisons$proxy_model, function(model_name) {
  fit <- models[[model_name]]
  covariance <- sandwich::vcovCL(fit, cluster = analysis$pwsid, type = "HC1")
  estimate <- coef(fit)[["potential_space_fraction_sampled"]] * space_sd
  se <- sqrt(covariance["potential_space_fraction_sampled", "potential_space_fraction_sampled"]) * space_sd
  data.frame(
    model = model_name, predictor_sd = space_sd, estimate_log_odds_per_sd = estimate,
    cluster_robust_se = se, odds_ratio_per_sd = exp(estimate),
    ci_low = exp(estimate - 1.96 * se), ci_high = exp(estimate + 1.96 * se),
    p_value = 2 * pnorm(-abs(estimate / se))
  )
}))

full <- models[["M5 + space + density"]]
beta <- coef(full)
covariance <- sandwich::vcovCL(full, cluster = analysis$pwsid, type = "HC1")
space_grid <- seq(quantile(analysis$potential_space_fraction_sampled, 0.05), quantile(analysis$potential_space_fraction_sampled, 0.95), length.out = 80)
adjusted <- bind_rows(lapply(space_grid, function(value) {
  new_data <- analysis
  new_data$potential_space_fraction_sampled <- value
  design <- model.matrix(delete.response(terms(full)), new_data)[, names(beta), drop = FALSE]
  response <- plogis(drop(design %*% beta))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  se <- sqrt(drop(t(gradient) %*% covariance %*% gradient))
  data.frame(
    potential_space_pct = value * 100, fitted_pct_green = mean(response) * 100,
    ci_low = 100 * max(0, mean(response) - 1.96 * se),
    ci_high = 100 * min(1, mean(response) + 1.96 * se)
  )
}))

output_dir <- here::here("outputs", "models", "landscape_space_experiment")
figure_dir <- file.path(dir_figures, "archived", "landscape_space_experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(components, file.path(output_dir, "overlap_proxy_components.csv"))
write_csv(proxy, file.path(output_dir, "supplier_region_proxy.csv"))
write_csv(candidates, file.path(output_dir, "sample_eligibility.csv"))
write_csv(analysis, file.path(output_dir, "model_sample.csv"))
write_csv(metrics, file.path(output_dir, "cv_metrics.csv"))
write_csv(cv_predictions, file.path(output_dir, "cv_predictions.csv"))
write_csv(comparisons, file.path(output_dir, "incremental_cv_comparison.csv"))
write_csv(coefficients, file.path(output_dir, "landscape_space_coefficients.csv"))
write_csv(adjusted, file.path(output_dir, "adjusted_space_predictions.csv"))
write_csv(count(analysis, region, name = "observations"), file.path(output_dir, "sample_by_region.csv"))
write_csv(
  data.frame(
    asset = "USGS/NLCD_RELEASES/2021_REL/NLCD/2021", native_resolution_m = 30,
    mask = "Landcover classes 21-24; impervious values 0-100 only",
    proxy = "1 - impervious_area_m2 / valid_developed_area_m2",
    aggregation = "Population-overlap-weighted block-group proxy over greenness-contributing overlaps",
    coverage_rule = "At least 99% valid NLCD developed-area coverage per block group and 99% proxy coverage of supplier green population weights",
    original_model_rows = nrow(retained), matched_model_rows = nrow(analysis),
    caveat = "Potential non-impervious space, not verified irrigable area; includes developed CII landscapes; NLCD 2021 versus imagery 2022",
    validation = "Exact retained PWS folds; six models use identical complete-case rows; retained model unchanged"
  ), file.path(output_dir, "analysis_definition.csv")
)
saveRDS(models, file.path(output_dir, "models.rds"))

labels <- metrics |>
  mutate(label = paste0("R2 = ", sprintf("%.3f", r_squared), "\nRMSE = ", sprintf("%.1f", rmse_percentage_points), " pp"))
comparison_plot <- ggplot(cv_predictions, aes(observed_pct, predicted_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_point(color = "#2c7fb8", alpha = 0.23, size = 1) +
  geom_text(data = labels, aes(x = 3, y = 97, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  coord_equal() +
  labs(
    x = "Observed developed-land greenness (%)", y = "Out-of-fold predicted greenness (%)",
    title = "Does potential landscape space improve prediction?",
    subtitle = paste(nrow(analysis), "matched records; <=132 GPCD; retained PWS folds"),
    caption = "Space is the non-impervious fraction within NLCD-developed pixels, not verified irrigable land.\nAll panels use identical records and folds; the retained primary model is unchanged."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "predicted_vs_observed.png"), comparison_plot, width = 9, height = 11, dpi = 180)

effect_plot <- ggplot() +
  geom_point(data = analysis, aes(potential_space_fraction_sampled * 100, pct_green_sampled), alpha = 0.15, color = "#777777", size = 1) +
  geom_ribbon(data = adjusted, aes(potential_space_pct, ymin = ci_low, ymax = ci_high), fill = "#2c7fb8", alpha = 0.18) +
  geom_line(data = adjusted, aes(potential_space_pct, fitted_pct_green), color = "#2c7fb8", linewidth = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    x = "Potential landscape-space fraction within developed land (%)", y = "Greenness within developed land (%)",
    title = "Adjusted greenness across potential landscape space",
    subtitle = "Income + space + density model, averaging over other observed covariates and regions",
    caption = "Curve shown over the 5th-95th percentiles; shading is a 95% PWS-clustered interval.\nNon-impervious space can be bare or unirrigated; associations are not causal."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "adjusted_greenness_vs_landscape_space.png"), effect_plot, width = 10, height = 7, dpi = 180)

cat("Original retained records:", nrow(retained), "matched complete records:", nrow(analysis), "\n")
print(metrics, width = Inf)
print(comparisons, row.names = FALSE)
print(coefficients, row.names = FALSE)
cat("Saved landscape-space experiment to", output_dir, "and", figure_dir, "\n")