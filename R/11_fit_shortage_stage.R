source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("sandwich", quietly = TRUE)) stop("Install the 'sandwich' package first.")
input <- file.path(dir_processed, "conservation", "stage_model_candidates_2022.csv")
if (!file.exists(input)) stop("Run R/05_download_conservation.R first.")
analysis <- read_csv(input, show_col_types = FALSE) |>
  filter(stage_eligible %in% TRUE, is.finite(peak_stage), between(outdoor_gpcd, 0, 132)) |>
  arrange(region, pwsid)
stopifnot(nrow(analysis) > 0L, !anyDuplicated(analysis[c("pwsid", "region")]))
stopifnot(all(analysis$peak_stage %in% 0:6), all(analysis$stage_known_months == 4L))
stage_values <- sort(unique(analysis$peak_stage))
stopifnot(0L %in% stage_values, length(stage_values) > 1L)
stage_labels <- paste("Stage", stage_values)
analysis <- analysis |>
  mutate(
    green_fraction = pct_green_sampled / 100,
    log_outdoor = log1p(outdoor_gpcd),
    log_income = log(income_proxy),
    region = factor(region),
    stage_category = factor(peak_stage, levels = stage_values, labels = stage_labels)
  )
stopifnot(all(is.finite(analysis$green_fraction)), all(between(analysis$green_fraction, 0, 1)))
continuous <- c("log_outdoor", "precip_gs_mm", "eto_gs_mm", "developed_fraction_sampled", "log_income")
stopifnot(all(vapply(analysis[continuous], function(values) all(is.finite(values)) && sd(values) > 0, logical(1))))

supplier_groups <- analysis |>
  group_by(stage_org_id) |>
  summarise(distinct_stages = n_distinct(peak_stage), peak_stage = first(peak_stage), .groups = "drop")
stopifnot(all(supplier_groups$distinct_stages == 1L))
stage_support <- analysis |>
  group_by(stage_category) |>
  summarise(
    observations = n(),
    pws_ids = n_distinct(pwsid),
    reporting_suppliers = n_distinct(stage_org_id),
    .groups = "drop"
  )
if (any(stage_support$reporting_suppliers < 2L)) {
  stop("A stage has fewer than two independent reporting suppliers; grouped validation cannot assess that category.")
}
region_support <- analysis |>
  group_by(region) |>
  summarise(reporting_suppliers = n_distinct(stage_org_id), .groups = "drop")
if (any(region_support$reporting_suppliers < 2L)) {
  stop("A region has fewer than two reporting suppliers; region-fixed-effect validation is not identifiable.")
}

formulas <- list(
  "M1: use + region" = green_fraction ~ log_outdoor + region,
  "M2: add climate" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + region,
  "M3: add development" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + region,
  "M4: add peak stage" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + stage_category + region,
  "M5: add income" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + stage_category + log_income + region,
  "Benchmark: income, no stage" = green_fraction ~ log_outdoor + precip_gs_mm + eto_gs_mm + developed_fraction_sampled + log_income + region
)
fit_model <- function(data, formula) {
  fit <- glm(formula, family = quasibinomial(link = "logit"), data = data)
  if (!fit$converged || any(!is.finite(coef(fit)))) stop("Model did not converge with a full-rank design.")
  fit
}
models <- lapply(formulas, function(formula) fit_model(analysis, formula))
covariances <- lapply(models, function(fit) sandwich::vcovCL(fit, cluster = analysis$stage_org_id, type = "HC1"))

set.seed(2022)
folds_valid <- FALSE
for (attempt in seq_len(1000)) {
  supplier_groups$fold <- NA_integer_
  for (stage in stage_values) {
    indices <- which(supplier_groups$peak_stage == stage)
    supplier_groups$fold[indices] <- sample(rep(seq_len(5), length.out = length(indices)))
  }
  analysis$fold <- supplier_groups$fold[match(analysis$stage_org_id, supplier_groups$stage_org_id)]
  folds_valid <- all(vapply(seq_len(5), function(fold) {
    train <- analysis[analysis$fold != fold, ]
    nrow(train) < nrow(analysis) && n_distinct(train$region) == n_distinct(analysis$region) &&
      n_distinct(train$stage_category) == length(stage_values) &&
      qr(model.matrix(formulas[["M5: add income"]], train))$rank == length(coef(models[["M5: add income"]]))
  }, logical(1)))
  if (folds_valid) break
}
if (!folds_valid) stop("Unable to construct five reporting-supplier folds with supported training categories.")

cv_predictions <- bind_rows(lapply(names(formulas), function(model_name) {
  predicted <- rep(NA_real_, nrow(analysis))
  for (fold in seq_len(5)) {
    test <- which(analysis$fold == fold)
    fit <- fit_model(analysis[analysis$fold != fold, ], formulas[[model_name]])
    predicted[test] <- 100 * predict(fit, newdata = analysis[test, ], type = "response")
  }
  data.frame(
    model = model_name, pwsid = analysis$pwsid, region = analysis$region,
    reporting_supplier = analysis$stage_org_id, fold = analysis$fold,
    observed_pct = analysis$pct_green_sampled, predicted_pct = predicted
  )
}))
stopifnot(all(is.finite(cv_predictions$predicted_pct)))
cv_predictions$model <- factor(cv_predictions$model, levels = names(formulas))
cv_metrics <- cv_predictions |>
  group_by(model) |>
  summarise(
    observations = n(), pws_ids = n_distinct(pwsid), reporting_suppliers = n_distinct(reporting_supplier),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    mae_percentage_points = mean(abs(observed_pct - predicted_pct)), .groups = "drop"
  )
comparison <- data.frame(
  reference_model = c("M3: add development", "Benchmark: income, no stage"),
  stage_model = c("M4: add peak stage", "M5: add income")
)
comparison$reference_r_squared <- cv_metrics$r_squared[match(comparison$reference_model, cv_metrics$model)]
comparison$stage_r_squared <- cv_metrics$r_squared[match(comparison$stage_model, cv_metrics$model)]
comparison$delta_r_squared <- comparison$stage_r_squared - comparison$reference_r_squared
comparison$reference_rmse <- cv_metrics$rmse_percentage_points[match(comparison$reference_model, cv_metrics$model)]
comparison$stage_rmse <- cv_metrics$rmse_percentage_points[match(comparison$stage_model, cv_metrics$model)]

stage_coefficients <- bind_rows(lapply(c("M4: add peak stage", "M5: add income"), function(model_name) {
  beta <- coef(models[[model_name]])
  terms <- grep("^stage_category", names(beta), value = TRUE)
  se <- sqrt(diag(covariances[[model_name]]))[terms]
  data.frame(
    model = model_name, stage = sub("^stage_category", "", terms),
    estimate_log_odds = unname(beta[terms]), cluster_robust_se = unname(se),
    odds_ratio_vs_stage_0 = exp(unname(beta[terms])),
    ci_low = exp(unname(beta[terms]) - 1.96 * unname(se)),
    ci_high = exp(unname(beta[terms]) + 1.96 * unname(se))
  )
}))
m5 <- models[["M5: add income"]]
beta <- coef(m5)
vcov_m5 <- covariances[["M5: add income"]]
stage_terms <- grep("^stage_category", names(beta), value = TRUE)
stage_cov <- vcov_m5[stage_terms, stage_terms, drop = FALSE]
wald <- drop(t(beta[stage_terms]) %*% solve(stage_cov, beta[stage_terms]))
joint_test <- data.frame(
  test = "M5 joint categorical peak-stage association",
  chi_squared = wald, degrees_of_freedom = length(stage_terms),
  p_value = pchisq(wald, df = length(stage_terms), lower.tail = FALSE)
)
adjusted <- bind_rows(lapply(stage_labels, function(stage_label) {
  new_data <- analysis
  new_data$stage_category <- factor(stage_label, levels = stage_labels)
  design <- model.matrix(delete.response(terms(m5)), new_data)
  design <- design[, names(beta), drop = FALSE]
  response <- plogis(drop(design %*% beta))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  se <- sqrt(drop(t(gradient) %*% vcov_m5 %*% gradient))
  data.frame(
    stage_category = stage_label, fitted_pct_green = 100 * mean(response),
    ci_low = 100 * max(0, mean(response) - 1.96 * se),
    ci_high = 100 * min(1, mean(response) + 1.96 * se)
  )
})) |>
  mutate(stage_category = factor(stage_category, levels = stage_labels))

output_dir <- here::here("outputs", "models", "shortage_stage_experiment")
figure_dir <- file.path(dir_figures, "archived", "shortage_stage_experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(analysis, file.path(output_dir, "model_sample.csv"))
write_csv(stage_support, file.path(output_dir, "stage_support.csv"))
write_csv(count(analysis, region, name = "observations"), file.path(output_dir, "sample_by_region.csv"))
write_csv(supplier_groups, file.path(output_dir, "reporting_supplier_folds.csv"))
write_csv(cv_predictions, file.path(output_dir, "cv_predictions.csv"))
write_csv(cv_metrics, file.path(output_dir, "cv_metrics.csv"))
write_csv(comparison, file.path(output_dir, "incremental_cv_comparison.csv"))
write_csv(stage_coefficients, file.path(output_dir, "stage_coefficients.csv"))
write_csv(joint_test, file.path(output_dir, "joint_stage_test.csv"))
write_csv(adjusted, file.path(output_dir, "adjusted_stage_predictions.csv"))
write_csv(
  data.frame(
    source_field = "DWR_STANDARD_LEVEL", documentation = "conservation_documentation.pdf, page 4",
    period = "2022-06-01 to 2022-09-30", seasonal_summary = "Maximum reported standardized stage; all four months required",
    predictor = "Categorical, stage 0 reference; not a numeric severity slope",
    validation = "5-fold CV grouped by ORG_ID and stratified by peak stage; supported training regions/categories required",
    fold_assignment_attempt = attempt,
    caveat = "Plan-response level, not verified irrigation restrictions, compliance, or causal policy effect. Sparse high stages."
  ), file.path(output_dir, "analysis_definition.csv")
)
saveRDS(models, file.path(output_dir, "models.rds"))

metric_labels <- cv_metrics |>
  mutate(label = paste0("R2 = ", sprintf("%.3f", r_squared), "\nRMSE = ", sprintf("%.1f", rmse_percentage_points), " pp"))
predicted_plot <- ggplot(cv_predictions, aes(observed_pct, predicted_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_point(color = "#2c7fb8", alpha = 0.3, size = 1.2) +
  geom_text(data = metric_labels, aes(x = 3, y = 97, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  coord_equal() +
  labs(
    x = "Observed developed-land greenness (%)", y = "Out-of-fold predicted greenness (%)",
    title = "Standardized shortage-stage experiment",
    subtitle = paste(nrow(analysis), "records; <=132 GPCD; reporting-supplier-grouped CV"),
    caption = "Peak June-September 2022 standardized stage is categorical. All panels use identical records and folds.\nStage indicates plan-response level, not policy effectiveness; the retained primary model is unchanged."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "predicted_vs_observed.png"), predicted_plot, width = 9, height = 11, dpi = 180)

coefficient_plot <- ggplot(stage_coefficients, aes(stage, odds_ratio_vs_stage_0, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#555555") +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  scale_y_log10() +
  coord_flip() +
  labs(
    x = NULL, y = "Odds ratio relative to stage 0 (log scale)", color = NULL,
    title = "Categorical peak-stage associations",
    subtitle = "Reporting-supplier-clustered 95% intervals; not causal policy effects",
    caption = "Stage 0 means no shortage level invoked, not absence of all conservation measures.\nHigh-stage categories have few reporting suppliers; see stage_support.csv."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "stage_coefficients.png"), coefficient_plot, width = 9, height = 6, dpi = 180)

stage_plot <- ggplot(adjusted, aes(stage_category, fitted_pct_green)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, color = "#2c7fb8") +
  geom_point(size = 2.5, color = "#2c7fb8") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Maximum reported stage, June-September 2022", y = "Adjusted greenness (%)",
    title = "Adjusted greenness by standardized shortage stage",
    subtitle = "M5 predictions averaged over observed covariates and regions",
    caption = "95% intervals clustered by reporting supplier. Stage categories are not equally spaced effects.\nSparse high stages and potential covariate extrapolation limit interpretation."
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "adjusted_greenness_by_stage.png"), stage_plot, width = 9, height = 6, dpi = 180)

print(stage_support, width = Inf)
print(cv_metrics, width = Inf)
print(comparison, row.names = FALSE)
cat("Joint stage test p =", format.pval(joint_test$p_value), "\nSaved experiment to", output_dir, "and", figure_dir, "\n")