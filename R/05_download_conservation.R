source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the 'readxl' package before processing conservation reports.")
}
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--refresh"))
output_dir <- file.path(dir_processed, "conservation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

resource_root <- paste0(
  "https://publish.data.ca.gov/dataset/",
  "c69ac02b-adfb-459a-bc58-bf69a8b572d2/resource/"
)
sources <- data.frame(
  file = c("conservation_actions_single.xlsx", "conservation_documentation.pdf", "conservation_monthly_supply_demand.csv"),
  url = c(
    paste0(resource_root, "70c0da06-4094-4110-8ac1-18d5c5883b21/download/cons_actions_single.xlsx"),
    paste0(resource_root, "54936a0d-09ce-46d3-a4e7-8808d42326c0/download/documentation-for-water-conservation-and-production-dataset_v1.1.pdf"),
    paste0(resource_root, "f4d50112-5fb5-4066-b45c-44696b10a49e/download/monthly_combined_dataset.csv")
  )
)
for (index in seq_len(nrow(sources))) {
  destination <- file.path(dir_raw, sources$file[index])
  if (!file.exists(destination) || "--refresh" %in% args) {
    temporary <- paste0(destination, ".tmp")
    download.file(sources$url[index], temporary, mode = "wb")
    if (grepl("xlsx$", destination)) {
      stopifnot("Single responses" %in% readxl::excel_sheets(temporary))
    } else if (grepl("pdf$", destination)) {
      stopifnot(rawToChar(readBin(temporary, "raw", n = 4)) == "%PDF")
    } else {
      headers <- names(read_csv(temporary, n_max = 0, show_col_types = FALSE))
      stopifnot(all(c("ORG_ID", "WATER_SYSTEM_ID", "DWR_STANDARD_LEVEL", "DWR_STANDARD_LEVEL_FLAG") %in% headers))
    }
    if (!file.copy(temporary, destination, overwrite = TRUE)) stop("Unable to save ", destination)
    unlink(temporary)
  }
}
sources$md5 <- unname(tools::md5sum(file.path(dir_raw, sources$file)))
sources$downloaded_utc <- format(file.info(file.path(dir_raw, sources$file))$mtime, tz = "UTC", usetz = TRUE)
sources$dataset_page <- "https://data.ca.gov/dataset/urws-conservation-supply-demand"
write_csv(sources, file.path(output_dir, "source_manifest.csv"))

raw <- readxl::read_excel(
  file.path(dir_raw, sources$file[1]),
  sheet = "Single responses",
  guess_max = 100000
)
required <- c(
  "Supplier Name", "Public Water System ID", "Reporting Month",
  "Enacted any measures", "Weekly watering restrictions",
  "Water Shortage Contingency Stage Invoked", "DWR Stage"
)
stopifnot(all(required %in% names(raw)))
stopifnot(inherits(raw[["Reporting Month"]], c("Date", "POSIXt")))
action_groups <- list(
  demand_reduction = c(
    "Raising rates", "Apply drought surcharges",
    "Reduced allocations (for agencies with budget-based rates)",
    "Residential water audits", "CII water audits", "Expanded existing rebate program",
    "Rationing", "Turf replacement/rebate", "Demand Reduction Other"
  ),
  use_restrictions = c(
    "Weekly watering restrictions", "Use type restrictions",
    "Excessive irrigation of outdoor landscapes",
    "Washing a motor vehicle with a hose not fitted with a shut-off nozzle",
    "Application of potable water directly to driveways or sidewalks",
    "Use of potable water in decorative water features",
    "The application of water to irrigate turf and ornamental landscapes during and within 48 hours after measurable rainfall",
    "Restrictions Other",
    "Not serving drinking water other than upon request in eating or drinking establishments",
    "Operators of hotels and motels providing guests with the option of choosing not to have towels and linens laundered daily",
    "Requiring water-efficient and drought-tolerant landscaping", "Industry Other"
  ),
  communication = c(
    "E-mails", "Paper mail", "Notifications via Customer App", "Website",
    "Articles/News releases", "Youtube", "Facebook", "Instagram", "Other Social Media",
    "Community events", "Door hanger", "Workshops", "Television", "Radio", "Billboard",
    "Paid Media Advertising", "Bus shelter", "Communication Other (Fill-in)"
  ),
  enforcement_actions = c(
    "Assigned a different rate tier", "Fine", "Penalty Other",
    "Notification via Door hanger", "Notification via Letter", "Notification via Other",
    "Notification via Phone call", "Notification via customer app"
  )
)
action_columns <- unlist(action_groups, use.names = FALSE)
stopifnot(length(action_columns) == 47L, !anyDuplicated(action_columns), all(action_columns %in% names(raw)))
action_flags <- lapply(raw[action_columns], function(values) toupper(trimws(as.character(values))))
stopifnot(all(vapply(action_flags, function(values) all(is.na(values) | values %in% c("", "Y")), logical(1))))
raw$checked_measure_count <- rowSums(as.data.frame(action_flags) == "Y", na.rm = TRUE)
write_csv(
  data.frame(action = action_columns, category = rep(names(action_groups), lengths(action_groups))),
  file.path(output_dir, "counted_actions.csv")
)
start_date <- as.Date("2022-06-01")
end_date <- as.Date("2022-10-01")
monthly <- raw |>
  transmute(
    reporting_supplier_key = trimws(as.character(.data[["Public Water System ID"]])),
    supplier_name = .data[["Supplier Name"]],
    reporting_date = as.Date(.data[["Reporting Month"]]),
    enacted_any_measures_raw = as.character(.data[["Enacted any measures"]]),
    checked_measure_count = .data$checked_measure_count,
    weekly_restriction_raw = trimws(as.character(.data[["Weekly watering restrictions"]])),
    local_stage_raw = as.character(.data[["Water Shortage Contingency Stage Invoked"]]),
    dwr_stage_raw = as.character(.data[["DWR Stage"]])
  ) |>
  filter(!is.na(reporting_date), reporting_date >= start_date, reporting_date < end_date) |>
  mutate(
    month = as.Date(format(reporting_date, "%Y-%m-01")),
    weekly_restriction_reported = if_else(toupper(weekly_restriction_raw) == "Y", 1, NA_real_),
    measure_response_status = case_when(
      toupper(trimws(enacted_any_measures_raw)) == "YES" & checked_measure_count > 0 ~ "reported_actions",
      toupper(trimws(enacted_any_measures_raw)) == "NO" & checked_measure_count == 0 ~ "explicit_no_actions",
      toupper(trimws(enacted_any_measures_raw)) == "NO" & checked_measure_count > 0 ~ "conflict_no_with_actions",
      TRUE ~ "unknown_or_empty_yes"
    ),
    reported_measure_count = if_else(
      measure_response_status %in% c("reported_actions", "explicit_no_actions"),
      as.numeric(checked_measure_count),
      NA_real_
    ),
    weekly_restriction_assumed_binary = case_when(
      !measure_response_status %in% c("reported_actions", "explicit_no_actions") ~ NA_real_,
      toupper(weekly_restriction_raw) == "Y" ~ 1,
      is.na(weekly_restriction_raw) | weekly_restriction_raw == "" ~ 0,
      TRUE ~ NA_real_
    ),
    reported_row = TRUE
  ) |>
  distinct()
source_reports <- monthly
write_csv(source_reports, file.path(output_dir, "supplier_reports_2022.csv"))
write_csv(count(source_reports, measure_response_status, name = "source_report_rows"),
          file.path(output_dir, "measure_response_audit.csv"))
monthly <- bind_rows(lapply(seq_len(nrow(source_reports)), function(index) {
  ids <- unique(trimws(strsplit(source_reports$reporting_supplier_key[index], ",", fixed = TRUE)[[1]]))
  stopifnot(all(grepl("^CA[0-9]{7}$", ids)))
  source_reports[rep(index, length(ids)), ] |>
    mutate(pwsid = ids)
}))
stopifnot(nrow(monthly) > 0L)
if (anyDuplicated(monthly[c("pwsid", "month")])) {
  stop("Duplicate supplier-month reports require review before aggregation.")
}
if (!all(is.na(monthly$weekly_restriction_raw) | toupper(monthly$weekly_restriction_raw) %in% c("", "Y"))) {
  stop("Unexpected restriction codes: review source definitions before assigning values.")
}
write_csv(monthly, file.path(output_dir, "monthly_reports_2022.csv"))

month_calendar <- data.frame(
  month = seq(start_date, by = "month", length.out = 4),
  days = c(30L, 31L, 31L, 30L)
)
grid <- expand.grid(
  pwsid = unique(monthly$pwsid),
  month = month_calendar$month,
  stringsAsFactors = FALSE
) |>
  left_join(month_calendar, by = "month") |>
  left_join(monthly, by = c("pwsid", "month"))
supplier <- grid |>
  group_by(pwsid) |>
  summarise(
    reporting_supplier_key = paste(sort(unique(na.omit(reporting_supplier_key))), collapse = ";"),
    reporting_supplier_keys = n_distinct(reporting_supplier_key, na.rm = TRUE),
    report_months = sum(reported_row %in% TRUE),
    positive_report_months = sum(weekly_restriction_reported == 1, na.rm = TRUE),
    unknown_restriction_months = sum(is.na(weekly_restriction_reported)),
    positive_report_days = sum(days[weekly_restriction_reported %in% 1]),
    unknown_restriction_days = sum(days[is.na(weekly_restriction_reported)]),
    known_measure_months = sum(is.finite(reported_measure_count)),
    unknown_measure_days = sum(days[!is.finite(reported_measure_count)]),
    measure_count_day_total = sum(days * reported_measure_count, na.rm = TRUE),
    assumed_binary_known_months = sum(!is.na(weekly_restriction_assumed_binary)),
    assumed_binary_yes_months = sum(weekly_restriction_assumed_binary == 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    binary_restriction_status = case_when(
      assumed_binary_known_months != 4L ~ "unknown_or_incomplete",
      assumed_binary_yes_months == 4L ~ "assumed_yes_all_months",
      assumed_binary_yes_months == 0L ~ "assumed_no_all_months",
      TRUE ~ "mixed_months"
    ),
    outdoor_restriction_assumed = case_when(
      binary_restriction_status == "assumed_yes_all_months" ~ 1,
      binary_restriction_status == "assumed_no_all_months" ~ 0,
      TRUE ~ NA_real_
    ),
    measure_months_complete = known_measure_months == 4L,
    mean_reported_measures = if_else(measure_months_complete, measure_count_day_total / 122, NA_real_),
    restriction_months_complete = unknown_restriction_months == 0L,
    weekly_restriction_fraction = if_else(
      restriction_months_complete,
      positive_report_days / 122,
      NA_real_
    )
  )
write_csv(supplier, file.path(output_dir, "supplier_summary_2022.csv"))

region_dirs <- sort(list.dirs(dir_processed, recursive = FALSE, full.names = TRUE))
region_dirs <- region_dirs[grepl("^hr_", basename(region_dirs))]
model_rows <- bind_rows(lapply(region_dirs, function(region_dir) {
  read_csv(file.path(region_dir, "pws_green_water_income.csv"), show_col_types = FALSE) |>
    mutate(region = basename(region_dir))
})) |>
  filter(
    water_data_complete %in% TRUE, is.finite(outdoor_gpcd), between(outdoor_gpcd, 0, 132),
    is.finite(pct_green_sampled), between(pct_green_sampled, 0, 100),
    is.finite(precip_gs_mm), precip_gs_mm >= 0,
    is.finite(eto_gs_mm), eto_gs_mm > 0,
    is.finite(developed_fraction_sampled), between(developed_fraction_sampled, 0, 1),
    is.finite(income_proxy), income_proxy > 0
  ) |>
  left_join(supplier, by = "pwsid") |>
  mutate(
    report_months = coalesce(report_months, 0L),
    positive_report_months = coalesce(positive_report_months, 0L),
    unknown_restriction_months = coalesce(unknown_restriction_months, 4L),
    restriction_months_complete = coalesce(restriction_months_complete, FALSE),
    known_measure_months = coalesce(known_measure_months, 0L),
    measure_months_complete = coalesce(measure_months_complete, FALSE),
    binary_restriction_status = coalesce(binary_restriction_status, "unknown_or_incomplete")
  )
write_csv(model_rows, file.path(output_dir, "model_candidates_2022.csv"))
coverage <- model_rows |>
  group_by(region) |>
  summarise(
    modeled_rows = n(),
    any_reporting_match = sum(report_months > 0),
    all_four_reports = sum(report_months == 4L),
    all_four_restriction_values_known = sum(restriction_months_complete),
    distinct_known_restriction_values = n_distinct(weekly_restriction_fraction, na.rm = TRUE),
    all_four_measure_counts_known = sum(measure_months_complete),
    distinct_known_measure_counts = n_distinct(mean_reported_measures, na.rm = TRUE),
    assumed_binary_yes = sum(outdoor_restriction_assumed == 1, na.rm = TRUE),
    assumed_binary_no = sum(outdoor_restriction_assumed == 0, na.rm = TRUE),
    assumed_binary_mixed = sum(binary_restriction_status == "mixed_months"),
    assumed_binary_unknown = sum(binary_restriction_status == "unknown_or_incomplete"),
    .groups = "drop"
  )
write_csv(coverage, file.path(output_dir, "coverage_by_region_2022.csv"))
known <- model_rows |>
  filter(measure_months_complete, is.finite(mean_reported_measures), reporting_supplier_keys == 1L)
model_ready <- n_distinct(known$mean_reported_measures) >= 2L
write_csv(
  data.frame(
    model_ready = model_ready,
    measure = "Day-weighted mean monthly count of reported conservation actions",
    known_model_rows = nrow(known),
    distinct_measure_values = n_distinct(known$mean_reported_measures),
    reason = if (model_ready) "Exploratory reported-activity count varies; not a restriction-intensity measure." else
      "Insufficient variation in complete four-month reported action counts."
  ),
  file.path(output_dir, "model_readiness.csv")
)
binary_candidates <- model_rows |>
  filter(is.finite(outdoor_restriction_assumed), reporting_supplier_keys == 1L)
write_csv(
  count(model_rows, region, binary_restriction_status, name = "supplier_region_records"),
  file.path(output_dir, "binary_restriction_coverage.csv")
)
write_csv(
  data.frame(
    model_ready = n_distinct(binary_candidates$outdoor_restriction_assumed) == 2L,
    eligible_model_rows = nrow(binary_candidates),
    assumed_yes = sum(binary_candidates$outdoor_restriction_assumed == 1),
    assumed_no = sum(binary_candidates$outdoor_restriction_assumed == 0),
    assumption = "Unchecked weekly restrictions mean No only in otherwise usable reports; constant four-month histories required."
  ),
  file.path(output_dir, "binary_restriction_readiness.csv")
)
cat("June-September source reports:", nrow(source_reports),
  "reporting suppliers:", n_distinct(source_reports$reporting_supplier_key),
  "expanded PWS IDs:", n_distinct(monthly$pwsid), "\n")
cat("Modeled suppliers with any report:", n_distinct(model_rows$pwsid[model_rows$report_months > 0]), "\n")
print(coverage, width = Inf)
print(count(source_reports, measure_response_status), width = Inf)
cat("Ready to model reported-action count:", model_ready, "\nSaved processing outputs to", output_dir, "\n")

stage_fields <- c(
  "ORG_ID", "SUPPLIER_NAME", "WATER_SYSTEM_ID", "REPORT_PERIOD_START_DATE",
  "REPORT_PERIOD_END_DATE", "DWR_STANDARD_LEVEL", "DWR_STANDARD_LEVEL_FLAG"
)
stage_raw <- read_csv(
  file.path(dir_raw, "conservation_monthly_supply_demand.csv"),
  col_select = all_of(stage_fields),
  col_types = cols(.default = col_character())
)
stage_source <- stage_raw |>
  transmute(
    stage_org_id = trimws(ORG_ID),
    stage_supplier_name = SUPPLIER_NAME,
    listed_pws_ids = trimws(WATER_SYSTEM_ID),
    month = as.Date(REPORT_PERIOD_START_DATE, format = "%Y-%m-%d"),
    period_end = as.Date(REPORT_PERIOD_END_DATE, format = "%Y-%m-%d"),
    standardized_stage_raw = trimws(DWR_STANDARD_LEVEL),
    standardized_stage_flag = trimws(DWR_STANDARD_LEVEL_FLAG)
  ) |>
  filter(!is.na(month), month >= start_date, month < end_date) |>
  mutate(
    stage_code_valid = !is.na(standardized_stage_raw) & grepl("^[0-6]( *\\([^)]*\\))?$", standardized_stage_raw),
    stage_valid = stage_code_valid & (is.na(standardized_stage_flag) | standardized_stage_flag == ""),
    stage = if_else(stage_valid, as.integer(substr(standardized_stage_raw, 1, 1)), NA_integer_),
    stage_report_present = TRUE
  ) |>
  distinct()
stopifnot(nrow(stage_source) > 0L, all(format(stage_source$month, "%d") == "01"))
stopifnot(all(!is.na(stage_source$stage_org_id) & nzchar(stage_source$stage_org_id)))
write_csv(stage_source, file.path(output_dir, "stage_source_reports_2022.csv"))
stage_monthly <- bind_rows(lapply(seq_len(nrow(stage_source)), function(index) {
  ids <- unique(trimws(strsplit(stage_source$listed_pws_ids[index], ",", fixed = TRUE)[[1]]))
  stopifnot(all(grepl("^CA[0-9]{7}$", ids)))
  stage_source[rep(index, length(ids)), ] |>
    mutate(pwsid = ids)
}))
if (anyDuplicated(stage_monthly[c("pwsid", "month")])) {
  stop("Standardized-stage data have duplicate PWS-month records; review source reporting units before aggregation.")
}
write_csv(stage_monthly, file.path(output_dir, "stage_monthly_pws_2022.csv"))
stage_grid <- expand.grid(
  pwsid = unique(stage_monthly$pwsid), month = month_calendar$month, stringsAsFactors = FALSE
) |>
  left_join(stage_monthly, by = c("pwsid", "month"))
stage_summary <- stage_grid |>
  group_by(pwsid) |>
  summarise(
    stage_report_months = sum(stage_report_present %in% TRUE),
    stage_known_months = sum(!is.na(stage)),
    stage_reporting_orgs = n_distinct(stage_org_id, na.rm = TRUE),
    stage_org_id = paste(sort(unique(na.omit(stage_org_id))), collapse = ";"),
    peak_stage = if (sum(!is.na(stage)) == 4L) max(stage) else NA_integer_,
    stage_changed_within_season = if (sum(!is.na(stage)) == 4L) n_distinct(stage) > 1L else NA,
    .groups = "drop"
  ) |>
  mutate(stage_eligible = stage_known_months == 4L & stage_reporting_orgs == 1L & !is.na(peak_stage))
write_csv(stage_summary, file.path(output_dir, "stage_supplier_summary_2022.csv"))
stage_candidates <- model_rows |>
  left_join(stage_summary, by = "pwsid") |>
  mutate(stage_eligible = coalesce(stage_eligible, FALSE))
write_csv(stage_candidates, file.path(output_dir, "stage_model_candidates_2022.csv"))
stage_coverage <- stage_candidates |>
  group_by(region) |>
  summarise(
    original_model_rows = n(),
    any_stage_report = sum(stage_report_months > 0, na.rm = TRUE),
    four_valid_stage_months = sum(stage_eligible),
    .groups = "drop"
  )
write_csv(stage_coverage, file.path(output_dir, "stage_coverage_by_region_2022.csv"))
stage_distribution <- stage_candidates |>
  filter(stage_eligible) |>
  group_by(peak_stage) |>
  summarise(
    observations = n(),
    pws_ids = n_distinct(pwsid),
    reporting_suppliers = n_distinct(stage_org_id),
    .groups = "drop"
  )
write_csv(stage_distribution, file.path(output_dir, "stage_peak_distribution_2022.csv"))
cat("Documented standardized-stage source reports:", nrow(stage_source),
    "reporting suppliers:", n_distinct(stage_source$stage_org_id), "\n")
print(stage_distribution, width = Inf)