source(here::here("R", "00_config.R"))
library(httr2)
library(sf)
library(dplyr)

regions_url <- paste0(
  "https://utility.arcgis.com/usrsvcs/servers/",
  "3ce65cc69e0048169396b3610d6256b0/rest/services/",
  "Boundaries/i03_Hydrologic_Regions/MapServer/0/query"
)

resp <- request(regions_url) |>
  req_url_query(
    where = "1=1", outFields = "HR_NAME", outSR = 4326,
    returnGeometry = "true", f = "geojson"
  ) |>
  req_timeout(120) |>
  req_perform()

regions <- read_sf(resp_body_string(resp)) |>
  group_by(HR_NAME) |>
  summarise(.groups = "drop") |>
  arrange(HR_NAME)

stopifnot(nrow(regions) == 10L, !anyNA(regions$HR_NAME))
write_sf(regions, file.path(dir_raw, "hydrologic_regions.gpkg"), delete_dsn = TRUE)
cat("Saved", nrow(regions), "hydrologic regions to", file.path(dir_raw, "hydrologic_regions.gpkg"), "\n")
print(regions$HR_NAME)