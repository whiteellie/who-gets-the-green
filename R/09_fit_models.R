source(here::here("R", "00_config.R"))
source(here::here("R", "landscape_space_helpers.R"))
library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Install the 'sandwich' package before fitting the clustered models.")
}
args <- commandArgs(trailingOnly = TRUE)
if (!all(args %in% c("--conservation", "--baseline", "--binary-restrictions"))) {
  stop("Default/--baseline fits the retained model; --conservation and --binary-restrictions reproduce archived experiments.")
}
if (length(unique(args)) > 1L) stop("Choose only one analysis mode.")
use_binary <- "--binary-restrictions" %in% args
use_conservation <- any(c("--conservation", "--binary-restrictions") %in% args)
conservation_term <- if (use_binary) "restriction_binary" else "measures_z"
if (use_conservation) {
  message("Reproducing an unsuccessful conservation experiment; not the retained primary model.")
}

region_dirs <- list.dirs(dir_processed, recursive = FALSE, full.names = TRUE)
region_dirs <- region_dirs[grepl("^hr_", basename(region_dirs))]
region_dirs <- sort(region_dirs)
stopifnot(length(region_dirs) == 10L)

region_data <- lapply(region_dirs, function(region_dir) {
  scope <- basename(region_dir)
  region_name <- tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope)))
  read_csv(
    file.path(region_dir, "pws_green_water_income.csv"),
    show_col_types = FALSE
  ) |>
    mutate(region = region_name)
})

all_complete <- bind_rows(region_data) |>
  filter(
    water_data_complete %in% TRUE,
    is.finite(outdoor_gpcd),
    outdoor_gpcd >= 0,
    is.finite(pct_green_sampled),
    between(pct_green_sampled, 0, 100),
    is.finite(precip_gs_mm),
    precip_gs_mm >= 0,
    is.finite(eto_gs_mm),
    eto_gs_mm > 0,
    is.finite(developed_fraction_sampled),
    between(developed_fraction_sampled, 0, 1),
    is.finite(income_proxy),
    income_proxy > 0
  )

all_complete <- mutate(all_complete, analysis_group = pwsid)
if (use_conservation) {
  conservation_path <- file.path(dir_processed, "conservation", "supplier_summary_2022.csv")
  if (!file.exists(conservation_path)) stop("Run R/05_download_conservation.R first.")
  conservation <- read_csv(conservation_path, show_col_types = FALSE)
  required_conservation <- c(
    "pwsid", "report_months", "reporting_supplier_key", "reporting_supplier_keys",
    if (use_binary) c("binary_restriction_status", "outdoor_restriction_assumed") else
      c("measure_months_complete", "mean_reported_measures")
  )
  stopifnot(all(required_conservation %in% names(conservation)), !anyDuplicated(conservation$pwsid))
  all_complete <- all_complete |>
    left_join(conservation, by = "pwsid") |>
    filter(
      report_months == 4L,
      reporting_supplier_keys == 1L, !is.na(reporting_supplier_key), nzchar(reporting_supplier_key)
    ) |>
    mutate(analysis_group = reporting_supplier_key)
  if (use_binary) {
    all_complete <- all_complete |>
      filter(
        binary_restriction_status %in% c("assumed_yes_all_months", "assumed_no_all_months"),
        outdoor_restriction_assumed %in% c(0, 1)
      )
  } else {
    all_complete <- all_complete |>
      filter(measure_months_complete %in% TRUE, is.finite(mean_reported_measures), between(mean_reported_measures, 0, 47))
  }
}
cv_group_ids <- sort(unique(all_complete$analysis_group))
gpcd_cutoff <- 132
if (!use_conservation) {
  space_data <- bind_rows(lapply(region_dirs, function(region_dir) {
    scope <- basename(region_dir)
    read_landscape_space(scope) |>
      mutate(region = tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope))))
  }))
  all_complete <- left_join(all_complete, space_data, by = c("pwsid", "region"))
  space_coverage <- select(all_complete, pwsid, region, outdoor_gpcd,
    potential_space_fraction_sampled, proxy_population_coverage, landscape_space_complete)
  stopifnot(all(all_complete$landscape_space_complete[all_complete$outdoor_gpcd <= gpcd_cutoff]))
  all_complete <- filter(all_complete, landscape_space_complete %in% TRUE)
}
primary_raw <- filter(all_complete, outdoor_gpcd <= gpcd_cutoff)
stopifnot(nrow(primary_raw) > 0L)
if (use_binary && n_distinct(primary_raw$outdoor_restriction_assumed) != 2L) {
  stop("Binary sensitivity requires both assumed Yes and No on the primary sample.")
}
if (use_conservation && !use_binary && n_distinct(primary_raw$mean_reported_measures) < 2L) {
  stop(paste(
    "Cannot estimate a reported-action association: no variation among complete four-month counts.",
    "Missing and conflicting monthly reports must remain unknown."
  ))
}
if (use_conservation && !use_binary) {
  measures_center <- mean(primary_raw$mean_reported_measures)
  measures_scale <- sd(primary_raw$mean_reported_measures)
}

gpcd_log_center <- mean(log1p(primary_raw$outdoor_gpcd))
gpcd_log_scale <- sd(log1p(primary_raw$outdoor_gpcd))
precip_center <- mean(primary_raw$precip_gs_mm)
precip_scale <- sd(primary_raw$precip_gs_mm)
eto_center <- mean(primary_raw$eto_gs_mm)
eto_scale <- sd(primary_raw$eto_gs_mm)
developed_center <- mean(primary_raw$developed_fraction_sampled)
developed_scale <- sd(primary_raw$developed_fraction_sampled)
income_log_center <- mean(log(primary_raw$income_proxy))
income_log_scale <- sd(log(primary_raw$income_proxy))
if (!use_conservation) {
  space_center <- mean(primary_raw$potential_space_fraction_sampled)
  space_scale <- sd(primary_raw$potential_space_fraction_sampled)
  stopifnot(is.finite(space_scale), space_scale > 0)
}
region_levels <- sort(unique(all_complete$region))

prepare_model_data <- function(data) {
  if (!use_conservation) {
    data <- mutate(data, space_z = (potential_space_fraction_sampled - space_center) / space_scale)
  }
  if (use_binary) {
    data <- mutate(data, restriction_binary = outdoor_restriction_assumed)
  } else if (use_conservation) {
    data <- mutate(data, measures_z = (mean_reported_measures - measures_center) / measures_scale)
  }
  data |>
    mutate(
      green_fraction = pct_green_sampled / 100,
      gpcd_log1p = log1p(outdoor_gpcd),
      income_log = log(income_proxy),
      gpcd_log_z = (gpcd_log1p - gpcd_log_center) / gpcd_log_scale,
      precip_z = (precip_gs_mm - precip_center) / precip_scale,
      eto_z = (eto_gs_mm - eto_center) / eto_scale,
      developed_z = (developed_fraction_sampled - developed_center) / developed_scale,
      income_log_z = (income_log - income_log_center) / income_log_scale,
      region = factor(region, levels = region_levels)
    )
}

analysis <- prepare_model_data(primary_raw)
unrestricted <- prepare_model_data(all_complete)

stopifnot(nrow(analysis) > 0L, n_distinct(analysis$region) == 10L)
stopifnot(all(analysis$green_fraction > 0 & analysis$green_fraction < 1))

fit_fractional_logit <- function(data, formula) {
  glm(formula, family = quasibinomial(link = "logit"), data = data)
}

full_formula <- green_fraction ~ gpcd_log_z + precip_z + eto_z + developed_z + income_log_z + region
full_model_name <- "M5: add income, <=132"
if (!use_conservation) full_formula <- update(full_formula, . ~ . + space_z)
if (use_conservation) full_formula <- update(full_formula, paste(". ~ . +", conservation_term))

models <- list(
  "M1: use + region, <=132" = fit_fractional_logit(
    analysis,
    green_fraction ~ gpcd_log_z + region
  ),
  "M2: add climate, <=132" = fit_fractional_logit(
    analysis,
    green_fraction ~ gpcd_log_z + precip_z + eto_z + region
  ),
  "M3: add development, <=132" = fit_fractional_logit(
    analysis,
    green_fraction ~ gpcd_log_z + precip_z + eto_z + developed_z + region
  )
)
if (use_conservation) {
  m4_name <- if (use_binary) "M4: assumed weekly restriction, <=132" else "M4: add reported actions, <=132"
  models[[m4_name]] <- fit_fractional_logit(
    analysis,
    update(full_formula, . ~ . - income_log_z)
  )
} else {
  models[["M4: add landscape space, <=132"]] <- fit_fractional_logit(
    analysis, update(full_formula, . ~ . - income_log_z)
  )
}
models[[full_model_name]] <- fit_fractional_logit(analysis, full_formula)
if (use_binary) {
  models[["Benchmark: income, no binary, <=132"]] <- fit_fractional_logit(
    analysis,
    update(full_formula, . ~ . - restriction_binary)
  )
}
sensitivity_name <- if (use_conservation) "Sensitivity: M5, all GPCD" else "Sensitivity: M5, proxy-covered"
models[[sensitivity_name]] <- fit_fractional_logit(unrestricted, full_formula)
model_data <- c(rep(list(analysis), length(models) - 1L), list(unrestricted))
names(model_data) <- names(models)
stopifnot(all(vapply(models, `[[`, logical(1), "converged")))
stopifnot(all(vapply(models, function(model) all(is.finite(coef(model))), logical(1))))
inference_unit <- if (use_conservation) "reporting supplier" else "PWS ID"

cluster_vcov <- function(model, data) {
  sandwich::vcovCL(model, cluster = data$analysis_group, type = "HC1")
}

predictor_terms <- c("gpcd_log_z", "precip_z", "eto_z", "developed_z", "income_log_z")
term_labels <- c(
  gpcd_log_z = "Outdoor use: log1p GPCD",
  precip_z = "Growing-season precipitation",
  eto_z = "Reference ET",
  developed_z = "Developed fraction",
  income_log_z = "Income proxy: log"
)
if (!use_conservation) {
  predictor_terms <- append(predictor_terms, "space_z", after = 4)
  term_labels <- append(term_labels, c(space_z = "Potential landscape-space fraction"), after = 4)
}
if (use_conservation) {
  predictor_terms <- append(predictor_terms, conservation_term, after = 4)
  conservation_label <- if (use_binary) "Weekly restriction: assumed Yes vs No" else "Reported conservation-action count"
  term_labels <- append(term_labels, setNames(conservation_label, conservation_term), after = 4)
}

coefficient_rows <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  data <- model_data[[model_name]]
  beta <- coef(model)
  robust_se <- sqrt(diag(cluster_vcov(model, data)))
  terms <- intersect(predictor_terms, names(beta))
  estimate <- unname(beta[terms])
  standard_error <- unname(robust_se[terms])
  data.frame(
    model = model_name,
    term = terms,
    estimate_log_odds = estimate,
    cluster_robust_se = standard_error,
    odds_ratio_per_sd = ifelse(terms == "restriction_binary", NA_real_, exp(estimate)),
    odds_ratio = exp(estimate),
    effect_unit = ifelse(terms == "restriction_binary", "Assumed Yes vs No", "One SD increase"),
    ci_low = exp(estimate - 1.96 * standard_error),
    ci_high = exp(estimate + 1.96 * standard_error),
    p_value = 2 * pnorm(-abs(estimate / standard_error))
  )
})
coefficients <- bind_rows(coefficient_rows) |>
  mutate(
    term_label = factor(unname(term_labels[term]), levels = unname(term_labels)),
    model = factor(model, levels = names(models))
  )

model_info <- bind_rows(lapply(names(models), function(model_name) {
  data <- model_data[[model_name]]
  data.frame(
    model = model_name,
    observations = nrow(data),
    unique_suppliers = n_distinct(data$pwsid),
    reporting_units = n_distinct(data$analysis_group),
    regions = n_distinct(data$region),
    max_outdoor_gpcd = max(data$outdoor_gpcd)
  )
}))

full_model <- models[[full_model_name]]
full_vcov <- cluster_vcov(full_model, analysis)
beta <- coef(full_model)
effect_contrast <- setNames(rep(0, length(beta)), names(beta))
effect_contrast[c("income_log_z", "gpcd_log_z")] <- 1
contrast_estimate <- sum(effect_contrast * beta)
contrast_se <- sqrt(drop(t(effect_contrast) %*% full_vcov %*% effect_contrast))
contrast_results <- data.frame(
  contrast = "signed income-plus-outdoor-use coefficient sum; not an absolute-magnitude comparison",
  signed_sum_log_odds_per_sd = contrast_estimate,
  cluster_robust_se = contrast_se,
  ci_low = contrast_estimate - 1.96 * contrast_se,
  ci_high = contrast_estimate + 1.96 * contrast_se,
  odds_ratio_ratio = exp(contrast_estimate),
  ratio_ci_low = exp(contrast_estimate - 1.96 * contrast_se),
  ratio_ci_high = exp(contrast_estimate + 1.96 * contrast_se),
  p_value_two_sided = 2 * pnorm(-abs(contrast_estimate / contrast_se))
)
prediction_max <- gpcd_cutoff
gpcd_values <- seq(0, prediction_max, length.out = 80)

prediction_rows <- lapply(gpcd_values, function(gpcd) {
  new_data <- analysis
  new_data$gpcd_log_z <- (log1p(gpcd) - gpcd_log_center) / gpcd_log_scale
  design <- model.matrix(delete.response(terms(full_model)), new_data)
  design <- design[, names(beta), drop = FALSE]
  response <- plogis(drop(design %*% beta))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  standard_error <- sqrt(as.numeric(t(gradient) %*% full_vcov %*% gradient))
  fitted_fraction <- mean(response)
  data.frame(
    outdoor_gpcd = gpcd,
    fitted_pct_green = 100 * fitted_fraction,
    ci_low = 100 * max(0, fitted_fraction - 1.96 * standard_error),
    ci_high = 100 * min(1, fitted_fraction + 1.96 * standard_error)
  )
})
predictions <- bind_rows(prediction_rows)

predict_marginal_at <- function(variable, value) {
  new_data <- analysis
  if (variable == "precip_gs_mm") {
    new_data$precip_z <- (value - precip_center) / precip_scale
  } else if (variable == "eto_gs_mm") {
    new_data$eto_z <- (value - eto_center) / eto_scale
  } else if (variable == "developed_fraction_sampled") {
    new_data$developed_z <- (value - developed_center) / developed_scale
  } else if (variable == "potential_space_fraction_sampled") {
    new_data$space_z <- (value - space_center) / space_scale
  } else if (variable == "income_proxy") {
    new_data$income_log_z <- (log(value) - income_log_center) / income_log_scale
  } else if (variable == "mean_reported_measures") {
    new_data$measures_z <- (value - measures_center) / measures_scale
  } else if (variable == "outdoor_restriction_assumed") {
    new_data$restriction_binary <- value
  }
  design <- model.matrix(delete.response(terms(full_model)), new_data)
  design <- design[, names(beta), drop = FALSE]
  response <- plogis(drop(design %*% beta))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  standard_error <- sqrt(as.numeric(t(gradient) %*% full_vcov %*% gradient))
  fitted_fraction <- mean(response)
  data.frame(
    fitted_pct_green = 100 * fitted_fraction,
    ci_low = 100 * max(0, fitted_fraction - 1.96 * standard_error),
    ci_high = 100 * min(1, fitted_fraction + 1.96 * standard_error)
  )
}

make_marginal_curve <- function(label, variable, values, plot_values = values) {
  bind_rows(lapply(seq_along(values), function(index) {
    predict_marginal_at(variable, values[index]) |>
      mutate(predictor = label, predictor_value = plot_values[index])
  }))
}

precip_grid <- seq(
  quantile(analysis$precip_gs_mm, 0.05),
  quantile(analysis$precip_gs_mm, 0.95),
  length.out = 60
)
eto_grid <- seq(
  quantile(analysis$eto_gs_mm, 0.05),
  quantile(analysis$eto_gs_mm, 0.95),
  length.out = 60
)
developed_grid <- seq(
  quantile(analysis$developed_fraction_sampled, 0.05),
  quantile(analysis$developed_fraction_sampled, 0.95),
  length.out = 60
)
income_grid <- exp(seq(
  quantile(log(analysis$income_proxy), 0.05),
  quantile(log(analysis$income_proxy), 0.95),
  length.out = 60
))

marginal_effects <- bind_rows(
  make_marginal_curve("Growing-season precipitation (mm)", "precip_gs_mm", precip_grid),
  make_marginal_curve("Reference ET (mm)", "eto_gs_mm", eto_grid),
  make_marginal_curve(
    "Developed fraction (%)",
    "developed_fraction_sampled",
    developed_grid,
    developed_grid * 100
  ),
  make_marginal_curve("Income proxy ($)", "income_proxy", income_grid)
)
marginal_order <- c(
  "Growing-season precipitation (mm)",
  "Reference ET (mm)",
  "Developed fraction (%)",
  "Income proxy ($)"
)
if (!use_conservation) {
  space_grid <- seq(quantile(analysis$potential_space_fraction_sampled, 0.05),
    quantile(analysis$potential_space_fraction_sampled, 0.95), length.out = 60)
  marginal_effects <- bind_rows(marginal_effects,
    make_marginal_curve("Potential landscape space (%)", "potential_space_fraction_sampled", space_grid, space_grid * 100))
  marginal_order <- append(marginal_order, "Potential landscape space (%)", after = 3)
}
if (use_conservation && !use_binary) {
  measures_grid <- seq(
    quantile(analysis$mean_reported_measures, 0.05),
    quantile(analysis$mean_reported_measures, 0.95),
    length.out = 60
  )
  marginal_effects <- bind_rows(
    marginal_effects,
    make_marginal_curve("Mean monthly reported actions", "mean_reported_measures", measures_grid)
  )
  marginal_order <- c(marginal_order, "Mean monthly reported actions")
}
marginal_effects <- mutate(marginal_effects, predictor = factor(predictor, levels = marginal_order))
marginal_points <- bind_rows(
  transmute(analysis, predictor = "Growing-season precipitation (mm)", predictor_value = precip_gs_mm, observed_pct = pct_green_sampled),
  transmute(analysis, predictor = "Reference ET (mm)", predictor_value = eto_gs_mm, observed_pct = pct_green_sampled),
  transmute(analysis, predictor = "Developed fraction (%)", predictor_value = developed_fraction_sampled * 100, observed_pct = pct_green_sampled),
  transmute(analysis, predictor = "Income proxy ($)", predictor_value = income_proxy, observed_pct = pct_green_sampled)
) |>
  mutate(predictor = factor(predictor, levels = marginal_order))

if (!use_conservation) {
  marginal_points <- bind_rows(marginal_points,
    transmute(analysis, predictor = "Potential landscape space (%)",
      predictor_value = potential_space_fraction_sampled * 100, observed_pct = pct_green_sampled)) |>
    mutate(predictor = factor(predictor, levels = marginal_order))
}
set.seed(2022)
if (use_conservation && !use_binary) {
  marginal_points <- bind_rows(
    marginal_points,
    transmute(analysis, predictor = "Mean monthly reported actions", predictor_value = mean_reported_measures, observed_pct = pct_green_sampled)
  ) |>
    mutate(predictor = factor(predictor, levels = marginal_order))
}
if (use_binary) {
  binary_adjusted_predictions <- bind_rows(lapply(c(0, 1), function(value) {
    predict_marginal_at("outdoor_restriction_assumed", value) |>
      mutate(restriction_binary = value, assumed_status = if (value == 1) "Reported Yes (Y)" else "Assumed No (unchecked)")
  }))
}
supplier_ids <- cv_group_ids
supplier_folds <- sample(rep(seq_len(5), length.out = length(supplier_ids)))
fold_map <- setNames(supplier_folds, supplier_ids)
cv_predictions <- bind_rows(lapply(names(models), function(model_name) {
  data <- model_data[[model_name]]
  data$cv_fold <- unname(fold_map[data$analysis_group])
  predicted_pct <- rep(NA_real_, nrow(data))
  model_formula <- formula(models[[model_name]])
  for (fold in seq_len(5)) {
    test_rows <- which(data$cv_fold == fold)
    training_data <- data[data$cv_fold != fold, ]
    test_data <- data[test_rows, ]
    cv_model <- fit_fractional_logit(training_data, model_formula)
    stopifnot(cv_model$converged)
    predicted_pct[test_rows] <- 100 * predict(cv_model, newdata = test_data, type = "response")
  }
  data.frame(
    model = model_name,
    pwsid = data$pwsid,
    analysis_group = data$analysis_group,
    region = data$region,
    cv_fold = data$cv_fold,
    observed_pct = data$pct_green_sampled,
    predicted_pct = predicted_pct
  )
}))
stopifnot(all(is.finite(cv_predictions$predicted_pct)))
cv_predictions$model <- factor(cv_predictions$model, levels = names(models))
cv_metrics <- cv_predictions |>
  group_by(model) |>
  summarise(
    method = paste("5-fold cross-validation grouped by", inference_unit),
    observations = n(),
    unique_suppliers = n_distinct(pwsid),
    reporting_units = n_distinct(analysis_group),
    folds = n_distinct(cv_fold),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    mae_percentage_points = mean(abs(observed_pct - predicted_pct)),
    .groups = "drop"
  ) |>
  mutate(
    model = factor(model, levels = names(models)),
    metric_label = paste0(
      "R2 = ", sprintf("%.3f", r_squared),
      "\nRMSE = ", sprintf("%.1f", rmse_percentage_points), " pp"
    )
  )
final_model_residuals <- NULL
if (!use_conservation) {
  final_model_residuals <- cv_predictions |>
    filter(as.character(model) == full_model_name) |>
    mutate(
      residual_pct = observed_pct - predicted_pct,
      cv_rmse_pct = sqrt(mean(residual_pct^2)),
      anomaly_cutoff_pct = 2 * cv_rmse_pct,
      anomaly_flag = abs(residual_pct) >= anomaly_cutoff_pct,
      residual_class = case_when(
        !anomaly_flag ~ "Within cutoff",
        residual_pct > 0 ~ "Greener than expected",
        TRUE ~ "Less green than expected"
      ),
      residual_class = factor(residual_class, levels = c(
        "Less green than expected", "Greener than expected", "Within cutoff"
      ))
    ) |>
    left_join(
      distinct(select(all_complete, pwsid, region, supplier_name)),
      by = c("pwsid", "region")
    ) |>
    arrange(desc(abs(residual_pct)))
  stopifnot(all(!is.na(final_model_residuals$supplier_name)))
}
progression_models <- names(models)[grepl("^M[1-5]:", names(models))]
incremental_cv_metrics <- cv_metrics |>
  filter(as.character(model) %in% progression_models) |>
  arrange(match(as.character(model), progression_models)) |>
  mutate(
    previous_model = lag(as.character(model)),
    incremental_r_squared = r_squared - lag(r_squared),
    incremental_r_squared_percentage_points = 100 * incremental_r_squared,
    added_predictor = case_when(
      grepl("^M2:", as.character(model)) ~ "Climate (precipitation + reference ET)",
      grepl("^M3:", as.character(model)) ~ "Developed fraction",
      grepl("^M4:", as.character(model)) & use_binary ~ "Assumed restrictions",
      grepl("^M4:", as.character(model)) & use_conservation ~ "Reported actions",
      grepl("^M4:", as.character(model)) ~ "Potential landscape space",
      grepl("^M5:", as.character(model)) ~ "Income",
      TRUE ~ "M1 baseline"
    )
  )
calibration <- cv_predictions |>
  group_by(model) |>
  mutate(prediction_decile = ntile(predicted_pct, 10)) |>
  group_by(model, prediction_decile) |>
  summarise(
    mean_predicted_pct = mean(predicted_pct),
    mean_observed_pct = mean(observed_pct),
    records = n(),
    .groups = "drop"
  ) |>
  arrange(model, prediction_decile)

model_dir <- if (use_binary) here::here("outputs", "models", "binary_restriction_sensitivity") else
  if (use_conservation) here::here("outputs", "models", "conservation_subset") else here::here("outputs", "models")
figure_root <- if (use_conservation) file.path(dir_figures, "archived") else dir_figures
figure_dir <- file.path(figure_root, if (use_binary) "binary_restriction_sensitivity_model" else
  if (use_conservation) "conservation_subset_model" else "all_regions_model")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
powerpoint_panel_dir <- file.path(dir_figures, "powerpoint_panels")
dir.create(powerpoint_panel_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(coefficients, file.path(model_dir, "all_regions_coefficients.csv"))
if (!use_conservation) write_csv(space_coverage, file.path(model_dir, "all_regions_landscape_space_coverage.csv"))
if (!use_conservation) write_csv(contrast_results, file.path(model_dir, "all_regions_m5_effect_contrast.csv"))
write_csv(model_info, file.path(model_dir, "all_regions_model_info.csv"))
write_csv(predictions, file.path(model_dir, "all_regions_adjusted_predictions.csv"))
write_csv(cv_predictions, file.path(model_dir, "all_regions_fractional_logit_cv_predictions.csv"))
write_csv(cv_metrics, file.path(model_dir, "all_regions_fractional_logit_cv_metrics.csv"))
if (!use_conservation) write_csv(final_model_residuals, file.path(model_dir, "all_regions_m5_cv_residuals.csv"))
write_csv(incremental_cv_metrics, file.path(model_dir, "incremental_cv_r_squared.csv"))
write_csv(calibration, file.path(model_dir, "all_regions_fractional_logit_calibration.csv"))
write_csv(marginal_effects, file.path(model_dir, "all_regions_m5_adjusted_marginal_effects.csv"))
write_csv(
  data.frame(
    variable = c("log1p_outdoor_gpcd", "precip_gs_mm", "eto_gs_mm", "developed_fraction_sampled",
      if (!use_conservation) "potential_space_fraction_sampled", "log_income_proxy"),
    mean = c(gpcd_log_center, precip_center, eto_center, developed_center, if (!use_conservation) space_center, income_log_center),
    sd = c(gpcd_log_scale, precip_scale, eto_scale, developed_scale, if (!use_conservation) space_scale, income_log_scale)
  ),
  file.path(model_dir, "all_regions_model_scaling.csv")
)
saveRDS(models, file.path(model_dir, "all_regions_fractional_logit_models.rds"))
if (use_conservation && !use_binary) {
  write_csv(
    data.frame(variable = "mean_reported_measures", mean = measures_center, sd = measures_scale),
    file.path(model_dir, "conservation_scaling.csv")
  )
}
full_model_description <- if (use_binary) "M5 (assumed binary restriction)" else
  if (use_conservation) "M5 (with reported actions)" else "M5 (landscape space + income)"
write_csv(count(analysis, region, name = "observations"), file.path(model_dir, "sample_by_region.csv"))
if (use_conservation) {
  sample_columns <- c(
    "pwsid", "region", "analysis_group", "outdoor_gpcd",
    if (use_binary) c("outdoor_restriction_assumed", "restriction_binary", "binary_restriction_status") else
      c("mean_reported_measures", "measures_z")
  )
  write_csv(
    select(analysis, all_of(sample_columns)),
    file.path(model_dir, "conservation_model_sample.csv")
  )
}
analysis_note <- if (use_binary) paste(
  "SENSITIVITY: unchecked weekly restrictions assumed No in usable reports.",
  "Mixed/incomplete histories excluded; this is not verified restriction status."
) else ""
if (use_binary) {
  write_csv(count(analysis, region, restriction_binary, name = "observations"), file.path(model_dir, "binary_group_counts.csv"))
  write_csv(data.frame(assumption = analysis_note, binary_predictor_standardized = FALSE), file.path(model_dir, "analysis_assumption.csv"))
  write_csv(binary_adjusted_predictions, file.path(model_dir, "binary_adjusted_predictions.csv"))
  comparison_pairs <- data.frame(
    reference_model = c("M3: add development, <=132", "Benchmark: income, no binary, <=132"),
    binary_model = c("M4: assumed weekly restriction, <=132", full_model_name)
  )
  comparison_pairs$reference_r_squared <- cv_metrics$r_squared[match(comparison_pairs$reference_model, cv_metrics$model)]
  comparison_pairs$binary_r_squared <- cv_metrics$r_squared[match(comparison_pairs$binary_model, cv_metrics$model)]
  comparison_pairs$delta_r_squared <- comparison_pairs$binary_r_squared - comparison_pairs$reference_r_squared
  comparison_pairs$reference_rmse <- cv_metrics$rmse_percentage_points[match(comparison_pairs$reference_model, cv_metrics$model)]
  comparison_pairs$binary_rmse <- cv_metrics$rmse_percentage_points[match(comparison_pairs$binary_model, cv_metrics$model)]
  write_csv(comparison_pairs, file.path(model_dir, "incremental_cv_comparison.csv"))
}

coefficient_plot <- ggplot(coefficients, aes(term_label, odds_ratio, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#555555") +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.65), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  coord_flip() +
  labs(
    x = NULL,
    y = if (use_binary) "Odds ratio: continuous per SD; binary assumed Yes vs No" else "Odds ratio per 1 SD increase (log scale)",
    color = NULL,
    title = if (use_binary) "Binary-restriction sensitivity: adjusted associations" else "Adjusted associations with developed-land greenness",
    subtitle = paste(nrow(analysis), "supplier-region records; uncertainty clustered by", inference_unit),
    caption = paste(
      if (use_binary) "M4 adds unverified assumed binary restrictions; M5 adds income. Continuous predictors are standardized; the binary is not." else
      if (use_conservation) "M1: water use + region; M2: + climate; M3: + development; M4: + reported actions; M5: + income." else
        "M1: use + region; M2: + climate; M3: + development; M4: + landscape space; M5: + income.",
      if (use_binary) paste(strwrap(analysis_note, width = 115), collapse = "\n") else
      if (use_conservation) "M4 adds reported-action count, not restriction intensity; M5 adds income. Unrestricted M5 is sensitivity." else
        "Space = non-impervious fraction of developed land, not irrigable area. Sensitivity: proxy-covered, no GPCD cutoff.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2))

ggsave(
  file.path(figure_dir, "model_coefficients.png"),
  coefficient_plot,
  width = 11,
  height = 7,
  dpi = 180
)

if (!use_conservation) {
  final_model_anomaly_plot <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#53645b", linewidth = 0.6) +
    geom_point(
      data = filter(final_model_residuals, !anomaly_flag),
      aes(predicted_pct, observed_pct),
      color = "#aeb9b1", alpha = 0.55, size = 1.6
    ) +
    geom_point(
      data = filter(final_model_residuals, anomaly_flag),
      aes(predicted_pct, observed_pct, fill = residual_class),
      shape = 21, color = "#1c2b25", stroke = 0.35, size = 2.8
    ) +
    scale_fill_manual(
      values = c("Less green than expected" = "#c45b42", "Greener than expected" = "#258778"),
      breaks = c("Greener than expected", "Less green than expected"),
      labels = c("Greener", "Less green"),
      name = NULL
    ) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    coord_equal() +
    labs(
      x = "Predicted greenness (%)",
      y = "Observed greenness (%)",
      title = "Observed vs predicted: final model",
      subtitle = paste0(
        nrow(final_model_residuals), " records | ", sum(final_model_residuals$anomaly_flag),
        " flagged at |residual| >= 2 x RMSE (", sprintf("%.1f", unique(final_model_residuals$anomaly_cutoff_pct)), " pp)"
      ),
      caption = paste(
        "Residual = observed - predicted. Predictions are held out by PWS; flagged points are a review screen, not a formal outlier test.",
        "Record IDs and residuals are in outputs/models/all_regions_m5_cv_residuals.csv.",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(figure_dir, "final_model_observed_vs_predicted_anomalies.png"),
    final_model_anomaly_plot,
    width = 10,
    height = 8,
    dpi = 180
  )
  residual_powerpoint_plot <- final_model_anomaly_plot +
    labs(title = NULL, subtitle = NULL, caption = NULL, fill = NULL) +
    guides(fill = "none") +
    theme_minimal(base_size = 18, base_family = "sans") +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "#66736a", linewidth = 0.5),
      axis.ticks = element_line(color = "#66736a", linewidth = 0.4),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(8, 2, 8, 8)
    )
  residual_legend_data <- data.frame(
    x = 0.12,
    y = c(0.65, 0.35),
    label = c("Greener", "Less green"),
    residual_class = factor(
      c("Greener than expected", "Less green than expected"),
      levels = levels(final_model_residuals$residual_class)
    )
  )
  residual_legend <- ggplot(residual_legend_data) +
    geom_point(aes(x, y, fill = residual_class), shape = 21, color = "#1c2b25", size = 4) +
    geom_text(aes(x = 0.3, y = y, label = label), hjust = 0, size = 18 / ggplot2::.pt, color = "#1c2b25") +
    scale_fill_manual(
      values = c("Less green than expected" = "#c45b42", "Greener than expected" = "#258778"),
      guide = "none"
    ) +
    scale_x_continuous(limits = c(0, 1.5), breaks = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(8, 8, 8, 0))
  standalone_scatter_panel <- patchwork::wrap_plots(
    residual_powerpoint_plot, residual_legend, ncol = 2, widths = c(0.74, 0.26)
  )
  ggsave(
    file.path(powerpoint_panel_dir, "03_final_model_observed_vs_predicted.png"),
    standalone_scatter_panel,
    width = 6.5,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

incremental_cv_metrics <- incremental_cv_metrics |>
  mutate(
    model_label = case_when(
      grepl("^M1:", as.character(model)) ~ "M1\nOutdoor use",
      grepl("^M2:", as.character(model)) ~ "M2\n+ Climate",
      grepl("^M3:", as.character(model)) ~ "M3\n+ Development",
      grepl("^M4:", as.character(model)) ~ "M4\n+ Landscape space",
      grepl("^M5:", as.character(model)) ~ "M5\n+ Income"
    ),
    model_label = factor(model_label, levels = unique(model_label)),
    bar_label = if_else(is.na(incremental_r_squared),
      sprintf("%.3f", r_squared),
      sprintf("%.3f\n(%+.3f)", r_squared, incremental_r_squared))
  )
cumulative_r_squared_plot <- ggplot(incremental_cv_metrics, aes(r_squared, model_label, fill = model_label)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = bar_label), hjust = -0.08, size = 18 / ggplot2::.pt, show.legend = FALSE) +
  scale_fill_manual(values = c(
    "M1\nOutdoor use" = "#c5df9b",
    "M2\n+ Climate" = "#9ccc75",
    "M3\n+ Development" = "#70b35a",
    "M4\n+ Landscape space" = "#3f9145",
    "M5\n+ Income" = "#176b3a"
  ), guide = "none") +
  scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1.2, 0.2),
    labels = scales::label_number(accuracy = 0.1), expand = expansion(mult = c(0, 0.01))) +
  scale_y_discrete(limits = rev(levels(incremental_cv_metrics$model_label)), drop = FALSE) +
  labs(
    x = "Cross-validated R²",
    y = NULL,
    title = "Cumulative cross-validated R² by model step",
    subtitle = paste(nrow(analysis), "records; M1 baseline followed by successive predictor additions"),
    caption = "Bar labels show cumulative R² and gain from the preceding step; same sample and PWS-grouped folds."
  ) +
  theme_minimal(base_size = 18) +
  theme(axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18))

ggsave(
  file.path(figure_dir, "cumulative_cv_r_squared.png"),
  cumulative_r_squared_plot,
  width = 9,
  height = 5.5,
  dpi = 180
)
if (!use_conservation) {
  ggsave(
    file.path(powerpoint_panel_dir, "04_cumulative_cv_r_squared.png"),
    cumulative_r_squared_plot +
      labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme_minimal(base_size = 18, base_family = "sans") +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "#66736a", linewidth = 0.5),
        axis.ticks = element_line(color = "#66736a", linewidth = 0.4),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.margin = margin(8, 8, 8, 8)
      ),
    width = 6.5,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

display_points <- analysis
ticks <- c(0, 1, 10, 25, 50, 100, 132)
ticks <- ticks[ticks <= prediction_max]
prediction_plot <- ggplot() +
  geom_point(
    data = display_points,
    aes(log1p(outdoor_gpcd), pct_green_sampled),
    color = "#777777",
    alpha = 0.28,
    size = 1.4
  ) +
  geom_ribbon(
    data = predictions,
    aes(log1p(outdoor_gpcd), ymin = ci_low, ymax = ci_high),
    fill = "#2c7fb8",
    alpha = 0.18
  ) +
  geom_line(
    data = predictions,
    aes(log1p(outdoor_gpcd), fitted_pct_green),
    color = "#2c7fb8",
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = log1p(ticks),
    labels = scales::label_comma()(ticks)
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Outdoor residential use (GPCD; log1p-scaled axis)",
    y = "Greenness within developed land (%)",
    title = "Adjusted greenness across outdoor water use",
    subtitle = paste0(full_model_description, ": marginal predictions averaged over observed covariates and regions"),
    caption = paste0(
      "Shading is a 95% interval clustered by ", inference_unit, ". The primary fit and plotted records include all ",
      "\n", nrow(analysis), " supplier-region records at or below ", gpcd_cutoff, " GPCD.",
      if (use_binary) paste0("\n", paste(strwrap(analysis_note, width = 105), collapse = "\n")) else ""
    )
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figure_dir, "adjusted_greenness_vs_outdoor_gpcd.png"),
  prediction_plot,
  width = 10,
  height = 7,
  dpi = 180
)

marginal_effects_plot <- ggplot() +
  geom_point(
    data = marginal_points,
    aes(predictor_value, observed_pct),
    color = "#777777",
    alpha = 0.12,
    size = 1
  ) +
  geom_ribbon(
    data = marginal_effects,
    aes(predictor_value, ymin = ci_low, ymax = ci_high),
    fill = "#2c7fb8",
    alpha = 0.18
  ) +
  geom_line(
    data = marginal_effects,
    aes(predictor_value, fitted_pct_green),
    color = "#2c7fb8",
    linewidth = 1
  ) +
  facet_wrap(~predictor, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Predictor value",
    y = "Greenness within developed land (%)",
    title = "Adjusted greenness across climate and context",
    subtitle = paste0(full_model_description, ": each curve varies one predictor over the observed data"),
    caption = paste0(
      "Other predictors and region are averaged over observed records; 95% intervals are clustered by ", inference_unit,
      ".\nAssociations are not causal.",
      if (use_binary) paste0("\n", paste(strwrap(analysis_note, width = 105), collapse = "\n")) else
        if (use_conservation) " Reported-action counts measure reporting breadth, not restriction effectiveness." else ""
    )
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figure_dir, "m5_adjusted_marginal_effects.png"),
  marginal_effects_plot,
  width = 11,
  height = if (use_conservation) 11 else 8,
  dpi = 180
)

predicted_observed_plot <- ggplot(cv_predictions, aes(observed_pct, predicted_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_point(color = "#2c7fb8", alpha = 0.32, size = 1.5) +
  geom_text(
    data = cv_metrics,
    aes(x = 3, y = 96, label = metric_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2
  ) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  coord_equal() +
  labs(
    x = "Observed greenness within developed land (%)",
    y = "Cross-validated predicted greenness (%)",
    title = "Observed vs predicted greenness by model",
    subtitle = paste("Five-fold CV grouped by", inference_unit, "| primary models use <=132 GPCD"),
    caption = paste(
      "Out-of-fold R² = 1 - SSE/SST; RMSE is in percentage points.",
      if (use_conservation) "All primary steps use the same sample and folds; sensitivity uses all complete GPCD records." else
        "Primary steps share samples and folds; sensitivity has no GPCD cutoff and requires >=99% proxy coverage.",
      if (use_binary) paste(strwrap(analysis_note, width = 90), collapse = "\n") else "",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figure_dir, "fractional_logit_models_predicted_vs_observed.png"),
  predicted_observed_plot,
  width = 8,
  height = if (use_binary) 11 else 8,
  dpi = 180
)

calibration_plot <- ggplot(calibration, aes(mean_predicted_pct, mean_observed_pct, group = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_line(color = "#2c7fb8", linewidth = 0.7) +
  geom_point(color = "#2c7fb8", size = 2) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  coord_equal() +
  labs(
    x = "Mean cross-validated predicted greenness (%)",
    y = "Mean observed greenness (%)",
    title = "Calibration by predicted decile",
    subtitle = paste("Out-of-fold predictions; five-fold CV grouped by", inference_unit),
    caption = paste(
      "Each point is an equal-count decile within a model.",
      "Points on the dashed line indicate agreement between mean predictions and observations.",
      if (use_binary) paste(strwrap(analysis_note, width = 90), collapse = "\n") else "",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figure_dir, "fractional_logit_models_calibration.png"),
  calibration_plot,
  width = 8,
  height = if (use_binary) 11 else 8,
  dpi = 180
)

if (use_binary) {
  binary_plot <- ggplot(binary_adjusted_predictions, aes(assumed_status, fitted_pct_green)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, color = "#2c7fb8") +
    geom_point(size = 3, color = "#2c7fb8") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      x = NULL,
      y = "Adjusted greenness within developed land (%)",
      title = "Binary restriction sensitivity: adjusted predictions",
      subtitle = "M5 predictions averaged over the same observed covariates and regions",
      caption = paste(
        paste(strwrap(analysis_note, width = 95), collapse = "\n"),
        "95% intervals are clustered by reporting supplier; associations are not causal.",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(figure_dir, "binary_adjusted_predictions.png"), binary_plot, width = 9, height = 6, dpi = 180)
}

cat("Primary supplier-region records (<=132 GPCD):", nrow(analysis), "\n")
if (!use_conservation) {
  cat(full_model_description, "signed income-plus-outdoor-use coefficient sum (not magnitudes):", round(contrast_estimate, 3),
    "95% CI", round(contrast_results$ci_low, 3), "to", round(contrast_results$ci_high, 3),
    "p =", format.pval(contrast_results$p_value_two_sided), "\n")
}
cat("Cross-validation grouped by", inference_unit, ":\n")
print(select(cv_metrics, model, observations, unique_suppliers, r_squared, rmse_percentage_points))
print(model_info, row.names = FALSE)
cat("Saved model outputs to", model_dir, "\nSaved model figures to", figure_dir, "\n")