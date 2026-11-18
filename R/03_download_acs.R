source(here::here("R", "00_config.R"))
library(tidycensus)
library(sf)
library(dplyr)

bg <- get_acs(
  geography = "block group", state = "CA", year = acs_year, survey = "acs5",
  variables = c(median_hh_income = "B19013_001", population = "B01003_001", housing_units = "B25001_001"),
  output = "wide", geometry = TRUE
) |>
  transmute(
    GEOID,
    median_hh_income = median_hh_incomeE,
    population = populationE,
    housing_units = housing_unitsE
  )

write_sf(bg, file.path(dir_raw, "acs_bg.gpkg"), delete_dsn = TRUE)
