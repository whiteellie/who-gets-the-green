library(here)

dir_raw <- here("data", "raw")
dir_interim <- here("data", "interim")
dir_processed <- here("data", "processed")
dir_figures <- here("outputs", "figures")
for (d in c(dir_raw, dir_interim, dir_processed, dir_figures)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# SWRCB drinking water system service area boundaries (ArcGIS layer); verify URL if it 404s
pws_boundary_url <- paste0(
  "https://gispublic.waterboards.ca.gov/portalserver/rest/services/",
  "Drinking_Water/California_Drinking_Water_System_Area_Boundaries/FeatureServer/0"
)

ckan_api <- "https://data.ca.gov/api/3/action"
water_use_query <- "urban water supplier monthly reporting residential water use"

acs_year <- 2023
# set once with tidycensus::census_api_key("KEY", install = TRUE)

year <- 2022 # latest year in the eAR monthly file (keep in sync with NDVI_YEAR in 04_greenness_gee.py)
indoor_baseline_months <- c(12, 1, 2) # winter months as indoor-use proxy

# Column map: fill in after inspecting the downloaded files
pws_id_col <- "SABL_PWSID" # boundaries
wu_pws_col <- "PWSID" # water use csv (eAR delivered water)
wu_date_col <- "DateStartOfMonth"
wu_res_gal_col <- "CALCULATED Quantity in Gallons"
wu_pop_col <- "Residential Population IN SDWIS"
wu_type_col <- "TypeName"
wu_res_types <- c("Single-family Residential", "Multi-family Residential")
