source(here::here("R", "00_config.R"))
library(httr2)
library(purrr)

# Step 1: list candidate datasets/resources on data.ca.gov
res <- request(paste0(ckan_api, "/package_search")) |>
  req_url_query(q = water_use_query, rows = 10) |>
  req_perform() |>
  resp_body_json()

walk(res$result$results, function(pkg) {
  cat("\n", pkg$name, "|", pkg$title, "\n")
  walk(pkg$resources, \(r) cat("   ", r$format, r$name, r$url, "\n"))
})

# Step 2: paste the chosen CSV resource URL here and run the rest
water_use_url <- NULL
if (!is.null(water_use_url)) {
  request(water_use_url) |>
    req_timeout(600) |>
    req_perform(path = file.path(dir_raw, "water_use_monthly.csv"))
}
