source(here::here("R", "00_config.R"))
source(here::here("R", "landscape_space_helpers.R"))
library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Install the 'sandwich' package before fitting the clustered models.")
}

args <- commandArgs(trailingOnly = TRUE)
scope <- if (length(args)) args[[1]] else "hr_san_francisco_bay"
if (length(args) > 1 || !grepl("^hr_[a-z0-9_]+$", scope)) {
  stop("Supply one hydrologic-region slug, e.g. hr_san_francisco_bay or hr_sacramento_river.")
}
region_name <- tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope)))
source_data <- read_csv(
  file.path(dir_processed, scope, "pws_green_water_income.csv"),
  show_col_types = FALSE
)
required <- c(
  "water_data_complete", "outdoor_gpcd", "pct_green_sampled", "income_proxy",
  "developed_fraction_sampled", "precip_gs_mm", "eto_gs_mm", "pwsid"
)
stopifnot(all(required %in% names(source_data)))

analysis_all <- source_data |>
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

gpcd_cutoff <- 132
cv_supplier_ids <- sort(unique(analysis_all$pwsid))
analysis_all <- left_join(analysis_all, read_landscape_space(scope), by = "pwsid")
space_coverage <- select(analysis_all, pwsid, outdoor_gpcd,
  potential_space_fraction_sampled, proxy_population_coverage, landscape_space_complete)
stopifnot(all(analysis_all$landscape_space_complete[analysis_all$outdoor_gpcd <= gpcd_cutoff]))
analysis_all <- filter(analysis_all, landscape_space_complete %in% TRUE)
primary_raw <- filter(analysis_all, outdoor_gpcd <= gpcd_cutoff)
stopifnot(nrow(primary_raw) > 0L, n_distinct(primary_raw$pwsid) == nrow(primary_raw))

scaling <- list(
  gpcd_log_center = mean(log1p(primary_raw$outdoor_gpcd)),
  gpcd_log_scale = sd(log1p(primary_raw$outdoor_gpcd)),
  precip_center = mean(primary_raw$precip_gs_mm),
  precip_scale = sd(primary_raw$precip_gs_mm),
  eto_center = mean(primary_raw$eto_gs_mm),
  eto_scale = sd(primary_raw$eto_gs_mm),
  developed_center = mean(primary_raw$developed_fraction_sampled),
  developed_scale = sd(primary_raw$developed_fraction_sampled),
  space_center = mean(primary_raw$potential_space_fraction_sampled),
  space_scale = sd(primary_raw$potential_space_fraction_sampled),
  income_log_center = mean(log(primary_raw$income_proxy)),
  income_log_scale = sd(log(primary_raw$income_proxy))
)
stopifnot(is.finite(scaling$space_scale), scaling$space_scale > 0)

prepare_data <- function(data) {
  data |>
    mutate(
      green_fraction = pct_green_sampled / 100,
      gpcd_log1p = log1p(outdoor_gpcd),
      income_log = log(income_proxy),
      gpcd_log_z = (gpcd_log1p - scaling$gpcd_log_center) / scaling$gpcd_log_scale,
      precip_z = (precip_gs_mm - scaling$precip_center) / scaling$precip_scale,
      eto_z = (eto_gs_mm - scaling$eto_center) / scaling$eto_scale,
      developed_z = (developed_fraction_sampled - scaling$developed_center) / scaling$developed_scale,
      space_z = (potential_space_fraction_sampled - scaling$space_center) / scaling$space_scale,
      income_log_z = (income_log - scaling$income_log_center) / scaling$income_log_scale
    )
}

analysis <- prepare_data(primary_raw)
unrestricted <- prepare_data(analysis_all)
stopifnot(all(analysis$green_fraction > 0 & analysis$green_fraction < 1))

fit_fractional_logit <- function(data, formula) {
  glm(formula, family = quasibinomial(link = "logit"), data = data)
}

formulas <- list(
  "M1: outdoor use" = green_fraction ~ gpcd_log_z,
  "M2: add climate" = green_fraction ~ gpcd_log_z + precip_z + eto_z,
  "M3: add development" = green_fraction ~ gpcd_log_z + precip_z + eto_z + developed_z,
  "M4: add landscape space" = green_fraction ~ gpcd_log_z + precip_z + eto_z + developed_z + space_z,
  "M5: add income" = green_fraction ~ gpcd_log_z + precip_z + eto_z + developed_z + space_z + income_log_z
)
models <- lapply(formulas, function(formula) fit_fractional_logit(analysis, formula))
models[["Sensitivity: M5, proxy-covered"]] <- fit_fractional_logit(unrestricted, formulas[["M5: add income"]])
model_formulas <- c(formulas, list("Sensitivity: M5, proxy-covered" = formulas[["M5: add income"]]))
model_data <- c(
  rep(list(analysis), length(formulas)),
  list(unrestricted)
)
names(model_data) <- names(models)
stopifnot(all(vapply(models, `[[`, logical(1), "converged")))
stopifnot(all(vapply(models, function(model) all(is.finite(coef(model))), logical(1))))

coefficient_rows <- lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  data <- model_data[[model_name]]
  beta <- coef(model)
  robust_vcov <- sandwich::vcovCL(model, cluster = data$pwsid, type = "HC1")
  robust_se <- sqrt(diag(robust_vcov))
  terms <- intersect(c("gpcd_log_z", "precip_z", "eto_z", "developed_z", "space_z", "income_log_z"), names(beta))
  estimates <- unname(beta[terms])
  standard_errors <- unname(robust_se[terms])
  data.frame(
    model = model_name,
    term = terms,
    estimate_log_odds = estimates,
    cluster_robust_se = standard_errors,
    odds_ratio_per_sd = exp(estimates),
    ci_low = exp(estimates - 1.96 * standard_errors),
    ci_high = exp(estimates + 1.96 * standard_errors),
    p_value = 2 * pnorm(-abs(estimates / standard_errors))
  )
})
term_labels <- c(
  gpcd_log_z = "Outdoor use: log1p GPCD",
  precip_z = "Growing-season precipitation",
  eto_z = "Reference ET",
  developed_z = "Developed fraction",
  space_z = "Potential landscape-space fraction",
  income_log_z = "Income proxy: log"
)
coefficients <- bind_rows(coefficient_rows) |>
  mutate(
    term_label = factor(unname(term_labels[term]), levels = unname(term_labels)),
    model = factor(model, levels = names(models))
  )

set.seed(2022)
supplier_ids <- cv_supplier_ids
supplier_folds <- sample(rep(seq_len(5), length.out = length(supplier_ids)))
fold_map <- setNames(supplier_folds, supplier_ids)
cv_predictions <- bind_rows(lapply(names(models), function(model_name) {
  data <- model_data[[model_name]]
  data$cv_fold <- unname(fold_map[data$pwsid])
  predicted_pct <- rep(NA_real_, nrow(data))
  for (fold in seq_len(5)) {
    test_rows <- which(data$cv_fold == fold)
    training_data <- data[data$cv_fold != fold, ]
    test_data <- data[test_rows, ]
    cv_model <- fit_fractional_logit(training_data, model_formulas[[model_name]])
    stopifnot(cv_model$converged)
    predicted_pct[test_rows] <- 100 * predict(cv_model, newdata = test_data, type = "response")
  }
  data.frame(
    model = model_name,
    pwsid = data$pwsid,
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
    observations = n(),
    unique_suppliers = n_distinct(pwsid),
    folds = n_distinct(cv_fold),
    r_squared = 1 - sum((observed_pct - predicted_pct)^2) / sum((observed_pct - mean(observed_pct))^2),
    rmse_percentage_points = sqrt(mean((observed_pct - predicted_pct)^2)),
    mae_percentage_points = mean(abs(observed_pct - predicted_pct)),
    .groups = "drop"
  )

m3 <- models[["M5: add income"]]
m3_vcov <- sandwich::vcovCL(m3, cluster = analysis$pwsid, type = "HC1")
marginal_outdoor <- lapply(seq(0, gpcd_cutoff, length.out = 80), function(gpcd) {
  new_data <- analysis
  new_data$gpcd_log_z <- (log1p(gpcd) - scaling$gpcd_log_center) / scaling$gpcd_log_scale
  design <- model.matrix(delete.response(terms(m3)), new_data)
  design <- design[, names(coef(m3)), drop = FALSE]
  response <- plogis(drop(design %*% coef(m3)))
  gradient <- colMeans(design * as.numeric(response * (1 - response)))
  standard_error <- sqrt(as.numeric(t(gradient) %*% m3_vcov %*% gradient))
  fitted_fraction <- mean(response)
  data.frame(
    outdoor_gpcd = gpcd,
    fitted_pct_green = 100 * fitted_fraction,
    ci_low = 100 * max(0, fitted_fraction - 1.96 * standard_error),
    ci_high = 100 * min(1, fitted_fraction + 1.96 * standard_error)
  )
}) |>
  bind_rows()

output_dir <- here::here("outputs", "models", scope)
figure_dir <- file.path(dir_figures, paste0(scope, "_model"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(coefficients, file.path(output_dir, "coefficients.csv"))
write_csv(space_coverage, file.path(output_dir, "landscape_space_coverage.csv"))
write_csv(cv_metrics, file.path(output_dir, "cv_metrics.csv"))
write_csv(cv_predictions, file.path(output_dir, "cv_predictions.csv"))
write_csv(marginal_outdoor, file.path(output_dir, "adjusted_outdoor_predictions.csv"))
write_csv(
  data.frame(
    variable = c("log1p_outdoor_gpcd", "precip_gs_mm", "eto_gs_mm", "developed_fraction_sampled", "potential_space_fraction_sampled", "log_income_proxy"),
    mean = unlist(scaling[c("gpcd_log_center", "precip_center", "eto_center", "developed_center", "space_center", "income_log_center")]),
    sd = unlist(scaling[c("gpcd_log_scale", "precip_scale", "eto_scale", "developed_scale", "space_scale", "income_log_scale")])
  ),
  file.path(output_dir, "scaling.csv")
)
saveRDS(models, file.path(output_dir, "models.rds"))

coefficient_plot <- ggplot(coefficients, aes(term_label, odds_ratio_per_sd, color = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#555555") +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.65), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  coord_flip() +
  labs(x = NULL, y = "Odds ratio per 1 SD increase (log scale)", color = NULL,
      title = paste(region_name, "adjusted associations"),
       subtitle = paste(nrow(analysis), "supplier records; uncertainty clustered by PWS")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(figure_dir, "model_coefficients.png"), coefficient_plot, width = 10, height = 6, dpi = 180)

prediction_plot <- ggplot(cv_predictions, aes(observed_pct, predicted_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#555555") +
  geom_point(color = "#2c7fb8", alpha = 0.35, size = 1.5) +
  geom_text(data = cv_metrics, aes(x = 3, y = 96, label = paste0("R2 = ", sprintf("%.3f", r_squared), "\nRMSE = ", sprintf("%.1f", rmse_percentage_points), " pp")), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~model, ncol = 2) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  coord_equal() +
  labs(x = "Observed greenness within developed land (%)", y = "Cross-validated predicted greenness (%)",
      title = paste(region_name, "observed vs predicted greenness"),
       subtitle = "Five-fold cross-validation grouped by PWS ID; primary models use <=132 GPCD",
       caption = paste(
         "M1: use; M2: + climate; M3: + development; M4: + landscape space; M5: + income.",
         "Sensitivity: no GPCD cutoff, >=99% proxy coverage. R2 = 1 - SSE/SST.",
         sep = "\n"
       )) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "predicted_vs_observed.png"), prediction_plot, width = 8, height = 8, dpi = 180)

outdoor_plot <- ggplot() +
  geom_point(data = analysis, aes(log1p(outdoor_gpcd), pct_green_sampled), color = "#777777", alpha = 0.3, size = 1.4) +
  geom_ribbon(data = marginal_outdoor, aes(log1p(outdoor_gpcd), ymin = ci_low, ymax = ci_high), fill = "#2c7fb8", alpha = 0.18) +
  geom_line(data = marginal_outdoor, aes(log1p(outdoor_gpcd), fitted_pct_green), color = "#2c7fb8", linewidth = 1) +
  scale_x_continuous(breaks = log1p(c(0, 1, 10, 25, 50, 100, 132)), labels = scales::label_comma()(c(0, 1, 10, 25, 50, 100, 132))) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Outdoor residential use (GPCD; log1p-scaled axis)", y = "Greenness within developed land (%)",
      title = paste(region_name, "adjusted greenness by outdoor use"),
      subtitle = "M5 (landscape space + income): marginal predictions over observed covariates",
       caption = paste("Shading: 95% PWS-clustered interval. Primary sample: <=132 GPCD.",
         "Landscape space is non-impervious developed land, not actual irrigable area.", sep = "\n")) +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "adjusted_greenness_vs_outdoor_gpcd.png"), outdoor_plot, width = 9, height = 6, dpi = 180)

cat(region_name, "complete records:", nrow(space_coverage), "\nProxy-covered sensitivity records:", nrow(analysis_all),
  "\nPrimary records <=132 GPCD:", nrow(analysis), "\n")
print(cv_metrics)
cat("Saved", region_name, "model outputs to", output_dir, "and figures to", figure_dir, "\n")