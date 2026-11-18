read_landscape_space <- function(scope, pws_ids = NULL, raw_dir = dir_raw, processed_dir = dir_processed) {
  stopifnot(grepl("^hr_[a-z0-9_]+$", scope))
  checkpoint <- file.path(raw_dir, paste0("landscape_space_bg_", scope, ".csv"))
  temporary <- sub("[.]csv$", ".tmp", checkpoint)
  if (!file.exists(checkpoint) && !file.exists(temporary)) {
    stop("Run python/06_landscape_space_gee.py for ", scope, " first.")
  }
  columns <- readr::cols(GEOID = "c", .default = readr::col_guess())
  nlcd <- if (file.exists(checkpoint)) readr::read_csv(checkpoint, col_types = columns) else NULL
  if (file.exists(temporary)) {
    recovered <- readr::read_csv(temporary, col_types = columns)
    if (!is.null(nlcd) && !all(nlcd$GEOID %in% recovered$GEOID)) {
      stop("Temporary landscape-space checkpoint is not a superset: ", temporary)
    }
    nlcd <- recovered
  }
  required <- c("GEOID", "potential_space_fraction", "impervious_coverage", "status")
  stopifnot(all(required %in% names(nlcd)), !anyDuplicated(nlcd$GEOID))
  overlap <- readr::read_csv(
    file.path(processed_dir, scope, "block_group_pws.csv"),
    col_types = columns,
    col_select = dplyr::all_of(c("GEOID", "pwsid", "pct_green_of_developed", "w"))
  )
  if (!is.null(pws_ids)) overlap <- dplyr::filter(overlap, pwsid %in% pws_ids)
  stopifnot(all(is.finite(overlap$w)), all(overlap$w > 0))
  overlap |>
    dplyr::left_join(nlcd, by = "GEOID") |>
    dplyr::mutate(
      green_valid = is.finite(pct_green_of_developed),
      space_valid = status %in% "ok" & is.finite(potential_space_fraction) &
        dplyr::between(potential_space_fraction, 0, 1) &
        is.finite(impervious_coverage) & impervious_coverage >= 0.99
    ) |>
    dplyr::group_by(pwsid) |>
    dplyr::summarise(
      proxy_population_coverage = if (sum(w[green_valid]) > 0)
        sum(w[green_valid & space_valid]) / sum(w[green_valid]) else NA_real_,
      potential_space_fraction_sampled = if (any(green_valid & space_valid))
        weighted.mean(potential_space_fraction[green_valid & space_valid], w[green_valid & space_valid]) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      landscape_space_complete = is.finite(proxy_population_coverage) & proxy_population_coverage >= 0.99 &
        is.finite(potential_space_fraction_sampled) &
        dplyr::between(potential_space_fraction_sampled, 0, 1)
    )
}