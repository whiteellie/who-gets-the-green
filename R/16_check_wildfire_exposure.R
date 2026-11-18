source(here::here("R", "00_config.R"))
library(sf)
library(dplyr)
library(readr)
library(httr2)
library(jsonlite)

fire_start_year <- 2015L
fire_end_year <- 2022L
search_distance_m <- 10000
page_size <- 1000L
fire_service <- paste0(
  "https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/arcgis/rest/services/",
  "California_Historic_Fire_Perimeters/FeatureServer/0/query"
)

residuals_all <- read_csv(
  file.path(here::here("outputs", "models"), "all_regions_m5_cv_residuals.csv"),
  show_col_types = FALSE
)
args <- commandArgs(trailingOnly = TRUE)
direction <- if (length(args)) args[[1]] else "less_green"
if (!direction %in% c("greener", "less_green")) {
  stop("Choose direction 'greener' or 'less_green'.")
}
residuals <- if (direction == "greener") {
  residuals_all |>
    filter(residual_pct > 0) |>
    arrange(desc(residual_pct)) |>
    slice_head(n = 3)
} else {
  residuals_all |>
    filter(residual_pct < 0) |>
    arrange(residual_pct) |>
    slice_head(n = 3)
}
residuals$direction <- direction
standout_ids <- unique(residuals$pwsid)
stopifnot(length(standout_ids) == 3L)

pws <- read_sf(file.path(dir_raw, "pws_boundaries.gpkg"), quiet = TRUE) |>
  filter(BOUNDARY_TYPE == "Water Service Area", SABL_PWSID %in% standout_ids) |>
  st_transform(3310) |>
  st_make_valid() |>
  transmute(pwsid = SABL_PWSID, supplier_name = WATER_SYSTEM_NAME) |>
  group_by(pwsid) |>
  summarise(supplier_name = first(supplier_name), do_union = TRUE, .groups = "drop")
stopifnot(setequal(standout_ids, pws$pwsid))

query_fire_page <- function(search_geometry, offset) {
  bbox <- st_bbox(st_transform(search_geometry, 3857))
  geometry <- toJSON(
    list(
      xmin = unname(bbox[["xmin"]]), ymin = unname(bbox[["ymin"]]),
      xmax = unname(bbox[["xmax"]]), ymax = unname(bbox[["ymax"]]),
      spatialReference = list(wkid = 3857)
    ),
    auto_unbox = TRUE, digits = 15
  )
  response <- request(fire_service) |>
    req_method("POST") |>
    req_body_form(
      where = paste0("YEAR_ >= ", fire_start_year, " AND YEAR_ <= ", fire_end_year, " AND STATE = 'CA'"),
      geometry = geometry,
      geometryType = "esriGeometryEnvelope",
      inSR = "3857",
      outSR = "3310",
      spatialRel = "esriSpatialRelIntersects",
      outFields = "YEAR_,FIRE_NAME,GIS_ACRES,ALARM_DATE,UNIT_ID",
      returnGeometry = "true",
      resultOffset = as.character(offset),
      resultRecordCount = as.character(page_size),
      f = "geojson"
    ) |>
    req_perform()
  response_text <- resp_body_string(response)
  response_object <- fromJSON(response_text, simplifyVector = FALSE)
  if (!is.null(response_object$error)) {
    stop("CAL FIRE perimeter query failed: ", response_object$error$message)
  }
  if (length(response_object$features) == 0L) return(NULL)
  page_path <- tempfile(fileext = ".geojson")
  writeLines(response_text, page_path, useBytes = TRUE)
  read_sf(page_path, quiet = TRUE) |>
    st_transform(3310)
}

query_nearby_fires <- function(system_row) {
  service_area <- st_geometry(system_row)
  search_area <- st_buffer(service_area, search_distance_m)
  pages <- list()
  offset <- 0L
  repeat {
    page <- query_fire_page(search_area, offset)
    if (is.null(page) || nrow(page) == 0L) break
    pages[[length(pages) + 1L]] <- page
    if (nrow(page) < page_size) break
    offset <- offset + page_size
  }
  if (length(pages) == 0L) return(NULL)
  fires <- bind_rows(pages) |>
    distinct(YEAR_, FIRE_NAME, GIS_ACRES, UNIT_ID, .keep_all = TRUE)
  distances_m <- as.numeric(st_distance(fires, service_area))
  intersects <- lengths(st_intersects(fires, service_area)) > 0L
  fires |>
    mutate(
      pwsid = system_row$pwsid,
      supplier_name = system_row$supplier_name,
      distance_to_service_area_km = distances_m / 1000,
      intersects_service_area = intersects
    ) |>
    filter(distance_to_service_area_km <= search_distance_m / 1000) |>
    select(
      pwsid, supplier_name, YEAR_, FIRE_NAME, GIS_ACRES, UNIT_ID,
      distance_to_service_area_km, intersects_service_area, geometry
    )
}

nearby_fires <- bind_rows(lapply(seq_len(nrow(pws)), function(index) {
  query_nearby_fires(pws[index, ])
}))
if (nrow(nearby_fires) == 0L) {
  stop("No fire perimeters were returned within the 10 km search envelopes.")
}

fire_details <- st_drop_geometry(nearby_fires) |>
  left_join(
    select(residuals, pwsid, region, residual_pct, observed_pct, predicted_pct),
    by = "pwsid"
  ) |>
  arrange(pwsid, distance_to_service_area_km)
fire_summary <- fire_details |>
  group_by(pwsid, supplier_name) |>
  summarise(
    regions = paste(sort(unique(region)), collapse = "; "),
    residual_pct = first(residual_pct),
    fire_perimeters_within_10km = n(),
    perimeters_intersecting_service_area = sum(intersects_service_area),
    nearest_perimeter_distance_km = min(distance_to_service_area_km),
    nearest_fire_names = paste(head(unique(FIRE_NAME), 5), collapse = "; "),
    fire_years = paste(sort(unique(YEAR_)), collapse = "; "),
    .groups = "drop"
  )
direct_fire_area <- bind_rows(lapply(seq_len(nrow(pws)), function(index) {
  service <- pws[index, ]
  direct <- filter(nearby_fires, pwsid == service$pwsid, intersects_service_area)
  service_area_m2 <- sum(as.numeric(st_area(service)))
  burned_area_m2 <- if (nrow(direct)) {
    intersection <- st_intersection(st_geometry(service), st_union(st_geometry(direct)))
    sum(as.numeric(st_area(intersection)))
  } else {
    0
  }
  data.frame(
    pwsid = service$pwsid,
    mapped_service_area_km2 = service_area_m2 / 1e6,
    unioned_direct_fire_area_km2 = burned_area_m2 / 1e6,
    direct_fire_area_pct = 100 * burned_area_m2 / service_area_m2
  )
}))
fire_summary <- left_join(fire_summary, direct_fire_area, by = "pwsid")

model_dir <- here::here("outputs", "models")
output_tag <- if (direction == "greener") "greener_than_expected" else "less_green_than_expected"
write_csv(fire_summary, file.path(model_dir, paste0(output_tag, "_wildfire_exposure_summary.csv")))
write_csv(fire_details, file.path(model_dir, paste0(output_tag, "_wildfire_perimeters_within_10km.csv")))
st_write(
  nearby_fires,
  file.path(model_dir, paste0(output_tag, "_wildfire_perimeters_within_10km.gpkg")),
  delete_dsn = TRUE,
  quiet = TRUE
)

cat("Direction:", direction, "| CAL FIRE historical perimeter exposure, years", fire_start_year, "to", fire_end_year,
  "; search distance", search_distance_m / 1000, "km\n")
print(fire_summary, n = Inf)