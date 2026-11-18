source(here::here("R", "00_config.R"))
library(httr2)
library(sf)

page <- 1000
offset <- 0
parts <- list()
repeat {
  resp <- request(paste0(pws_boundary_url, "/query")) |>
    req_url_query(
      where = "1=1", outFields = "*", outSR = 4326, f = "geojson",
      resultOffset = offset, resultRecordCount = page
    ) |>
    req_timeout(120) |>
    req_perform()
  batch <- read_sf(resp_body_string(resp))
  if (nrow(batch) == 0) break
  parts[[length(parts) + 1]] <- batch
  offset <- offset + nrow(batch)
}
pws <- do.call(rbind, parts)
print(names(pws))
write_sf(pws, file.path(dir_raw, "pws_boundaries.gpkg"), delete_dsn = TRUE)
