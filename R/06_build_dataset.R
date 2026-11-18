source(here::here("R", "00_config.R"))
library(sf)
library(dplyr)
library(readr)
library(lubridate)

args <- commandArgs(trailingOnly = TRUE)
scope <- if (length(args)) args[[1]] else "statewide"
region_scope <- grepl("^hr_[a-z0-9_]+$", scope)
county_scope <- grepl("^06[0-9]{3}$", scope)
stopifnot(scope == "statewide" || county_scope || region_scope)
green_file <- if (scope == "statewide") "greenness_bg.csv" else paste0("greenness_bg_", scope, ".csv")
output_dir <- file.path(dir_processed, scope)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

region <- NULL
if (region_scope) {
  regions <- read_sf(file.path(dir_raw, "hydrologic_regions.gpkg"))
  region_slug <- gsub("[^a-z0-9]+", "_", tolower(regions$HR_NAME))
  region <- regions[region_slug == sub("^hr_", "", scope), ] |>
    st_transform(3310) |>
    st_make_valid()
  stopifnot(nrow(region) == 1L)
}

wu <- read_csv(
  file.path(dir_raw, "water_use_monthly.csv"),
  col_select = all_of(c(wu_pws_col, wu_date_col, wu_type_col, wu_res_gal_col, wu_pop_col)),
  show_col_types = FALSE
) |>
  filter(.data[[wu_type_col]] %in% wu_res_types) |>
  group_by(pwsid = .data[[wu_pws_col]], date = mdy(.data[[wu_date_col]])) |>
  summarise(
    gal = if (all(is.na(.data[[wu_res_gal_col]]))) NA_real_ else sum(.data[[wu_res_gal_col]], na.rm = TRUE),
    pop = if (all(is.na(.data[[wu_pop_col]]))) NA_real_ else max(.data[[wu_pop_col]], na.rm = TRUE),
    reported_categories = sum(!is.na(.data[[wu_res_gal_col]])),
    rows = n(),
    .groups = "drop"
  ) |>
  filter(lubridate::year(date) == .env$year, is.finite(pop), pop > 0) |>
  mutate(
    days = days_in_month(date),
    gpcd = if_else(reported_categories == 2L & rows == 2L & gal >= 0, gal / days / pop, NA_real_)
  )

outdoor <- wu |>
  group_by(pwsid) |>
  summarise(
    valid_months = sum(is.finite(gpcd)),
    winter_months = sum(is.finite(gpcd) & month(date) %in% indoor_baseline_months),
    total_gpcd = weighted.mean(gpcd, days, na.rm = TRUE),
    indoor_gpcd = weighted.mean(gpcd[month(date) %in% indoor_baseline_months], days[month(date) %in% indoor_baseline_months], na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    water_data_complete = valid_months == 12L & winter_months == 3L & is.finite(total_gpcd) & total_gpcd > 0,
    water_quality = case_when(
      valid_months != 12L | winter_months != 3L ~ "incomplete_months",
      !is.finite(total_gpcd) | total_gpcd <= 0 ~ "zero_or_invalid_annual_use",
      TRUE ~ "usable_proxy"
    ),
    outdoor_gpcd = if_else(water_data_complete, pmax(total_gpcd - indoor_gpcd, 0), NA_real_)
  )

pws <- read_sf(file.path(dir_raw, "pws_boundaries.gpkg")) |>
  filter(BOUNDARY_TYPE == "Water Service Area", !is.na(.data[[pws_id_col]])) |>
  st_transform(3310) |>
  st_make_valid() |>
  select(pwsid = all_of(pws_id_col), supplier_name = WATER_SYSTEM_NAME) |>
  group_by(pwsid) |>
  summarise(supplier_name = first(supplier_name), .groups = "drop")
pws$service_area_m2 <- as.numeric(st_area(pws))

green <- read_csv(file.path(dir_raw, green_file), col_types = cols(GEOID = "c", .default = "d"))
stopifnot(!anyDuplicated(green$GEOID), nrow(green) > 0)
stopifnot(all(is.na(green$pct_green_of_developed) | between(green$pct_green_of_developed, 0, 1)))
if (county_scope) stopifnot(all(startsWith(green$GEOID, scope)))
if (region_scope) {
  climate <- read_csv(
    file.path(dir_raw, paste0("climate_bg_", scope, ".csv")),
    col_types = cols(GEOID = "c", .default = "d")
  )
  stopifnot(!anyDuplicated(climate$GEOID))
  stopifnot(all(c("precip_gs_mm", "eto_gs_mm") %in% names(climate)))
  stopifnot(all(green$GEOID %in% climate$GEOID))
  green <- left_join(green, climate, by = "GEOID")
  stopifnot(all(is.finite(green$precip_gs_mm)), all(is.finite(green$eto_gs_mm)))
}

bg <- read_sf(file.path(dir_raw, "acs_bg.gpkg")) |>
  semi_join(green, by = "GEOID") |>
  st_transform(3310) |>
  st_make_valid() |>
  inner_join(green, by = "GEOID")
bg$bg_area <- as.numeric(st_area(bg))
if (region_scope) bg <- st_intersection(bg, region |> select(HR_NAME))
pws <- st_filter(pws, bg)

overlap <- st_intersection(bg, pws)
overlap$overlap_area_m2 <- as.numeric(st_area(overlap))
overlap <- overlap |>
  mutate(w = overlap_area_m2 / bg_area * coalesce(population, 0)) |>
  st_drop_geometry() |>
  filter(w > 0)

wavg <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
}

pws_tbl <- overlap |>
  group_by(pwsid) |>
  summarise(
    supplier_name = first(supplier_name),
    income_proxy = wavg(median_hh_income, w),
    pct_green_sampled = wavg(pct_green_of_developed, w) * 100,
    developed_fraction_sampled = wavg(developed, w),
    precip_gs_mm = if (region_scope) wavg(precip_gs_mm, w) else NA_real_,
    eto_gs_mm = if (region_scope) wavg(eto_gs_mm, w) else NA_real_,
    sampled_population_est = sum(w),
    sampled_block_groups = n_distinct(GEOID),
    sampled_service_area_pct = 100 * sum(overlap_area_m2) / first(service_area_m2),
    .groups = "drop"
  )

out <- pws_tbl |>
  left_join(outdoor, by = "pwsid") |>
  mutate(
    water_year = year,
    acs_end_year = acs_year,
    green_per_outdoor_gpcd = if_else(outdoor_gpcd > 0, pct_green_sampled / outdoor_gpcd, NA_real_)
  )

block_group_columns <- c(
  "GEOID", "pwsid", "supplier_name", "median_hh_income", "population",
  "pct_green_of_developed", "developed", "overlap_area_m2", "bg_area", "w",
  if (region_scope) c("precip_gs_mm", "eto_gs_mm")
)
block_groups <- overlap |>
  select(all_of(block_group_columns)) |>
  left_join(outdoor, by = "pwsid")
write_csv(out, file.path(output_dir, "pws_green_water_income.csv"))
write_csv(block_groups, file.path(output_dir, "block_group_pws.csv"))
stopifnot(!anyDuplicated(out$pwsid), all(out$sampled_service_area_pct <= 100.01))
cat("Sampled block groups:", nrow(bg), "\nMatched suppliers:", nrow(out),
    "\nSuppliers with complete monthly water data:", sum(out$water_data_complete %in% TRUE), "\n")
print(out, width = Inf)
