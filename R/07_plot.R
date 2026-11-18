source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
scope <- if (length(args)) args[[1]] else "statewide"
max_outdoor_gpcd <- if (length(args) >= 2) as.numeric(args[[2]]) else Inf
region_scope <- grepl("^hr_[a-z0-9_]+$", scope)
county_scope <- grepl("^06[0-9]{3}$", scope)
stopifnot(scope == "statewide" || county_scope || region_scope)
stopifnot(length(max_outdoor_gpcd) == 1L, !is.na(max_outdoor_gpcd), max_outdoor_gpcd > 0)
figure_scope <- if (is.finite(max_outdoor_gpcd)) paste0(scope, "_gpcd_max_", max_outdoor_gpcd) else scope
output_dir <- file.path(dir_figures, figure_scope)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
df_all <- read_csv(file.path(dir_processed, scope, "pws_green_water_income.csv"), show_col_types = FALSE) |>
  filter(water_data_complete, is.finite(outdoor_gpcd), is.finite(pct_green_sampled), is.finite(income_proxy)) |>
  arrange(income_proxy)
excluded_high_gpcd <- sum(df_all$outdoor_gpcd > max_outdoor_gpcd)
df <- df_all |>
  filter(outdoor_gpcd <= max_outdoor_gpcd) |>
  mutate(
    label_y = pmin(98, pmax(2, pct_green_sampled + rep(c(3, -3), length.out = n()))),
    label_pwsid = if_else(min_rank(desc(outdoor_gpcd)) <= 10L, pwsid, NA_character_)
  )
stopifnot(nrow(df) > 0)
caption <- paste(
  if (region_scope) "Hydrologic-region analysis: sampled greenness vs supplier-wide water-use proxy; not an efficiency ranking." else "Pilot only: sampled greenness vs supplier-wide water-use proxy; not an efficiency ranking.",
  if (is.finite(max_outdoor_gpcd)) paste0(excluded_high_gpcd, " suppliers above ", max_outdoor_gpcd, " GPCD excluded from scatterplots; source data unchanged.") else NULL,
  "Water and imagery: 2022. ACS: 2019-2023. No climate adjustment.",
  "Income is an overlap-population-weighted mean of block-group medians, not a supplier median.",
  sep = "\n"
)
theme_set(theme_minimal(base_size = 11))

p1 <- ggplot(df, aes(outdoor_gpcd, pct_green_sampled, fill = pct_green_sampled, size = income_proxy)) +
  geom_point(shape = 21, color = "black", stroke = 0.45, alpha = 0.75) +
  geom_text(data = filter(df, !is.na(label_pwsid)), aes(y = label_y, label = label_pwsid),
            color = "#303030", size = 3, check_overlap = TRUE) +
  scale_fill_gradientn(
    colors = c("#d95f02", "#f7e225", "#1a9850"), limits = c(0, 100),
    labels = scales::label_percent(scale = 1), name = "Sampled greenness"
  ) +
  scale_size_continuous(name = "Income proxy", labels = scales::label_dollar(), range = c(1.5, 11)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Estimated outdoor residential use (gallons/person/day)", y = "Sampled developed-land greenness (%)",
    title = "Who gets the green?",
    subtitle = paste("Scope:", scope, "|", nrow(df), "suppliers with usable water records"), caption = caption
  )
ggsave(file.path(output_dir, "green_vs_outdoor_use.png"), p1, width = 10, height = 7, dpi = 180)

p2 <- ggplot(df, aes(income_proxy, pct_green_sampled, fill = pct_green_sampled, size = outdoor_gpcd)) +
  geom_point(shape = 21, color = "black", stroke = 0.45, alpha = 0.75) +
  geom_text(data = filter(df, !is.na(label_pwsid)), aes(y = label_y, label = label_pwsid),
            color = "#303030", size = 3, check_overlap = TRUE) +
  scale_fill_gradientn(
    colors = c("#d95f02", "#f7e225", "#1a9850"), limits = c(0, 100),
    labels = scales::label_percent(scale = 1), name = "Sampled greenness"
  ) +
  scale_size_continuous(name = "Outdoor GPCD", range = c(2, 7)) +
  scale_x_continuous(labels = scales::label_dollar()) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Sampled household income proxy", y = "Sampled developed-land greenness (%)",
      title = "Greenness and income across suppliers", caption = caption)
ggsave(file.path(output_dir, "green_vs_income.png"), p2, width = 10, height = 7, dpi = 180)

green_file <- if (scope == "statewide") "greenness_bg.csv" else paste0("greenness_bg_", scope, ".csv")
green <- read_csv(file.path(dir_raw, green_file), col_types = cols(GEOID = "c", .default = "d"))
geography <- sf::read_sf(file.path(dir_raw, "acs_bg.gpkg"))
if (county_scope) geography <- filter(geography, startsWith(GEOID, scope))
if (region_scope) {
  regions <- sf::read_sf(file.path(dir_raw, "hydrologic_regions.gpkg"))
  region_slug <- gsub("[^a-z0-9]+", "_", tolower(regions$HR_NAME))
  region <- regions[region_slug == sub("^hr_", "", scope), ]
  stopifnot(nrow(region) == 1L)
  geography <- sf::st_intersection(geography, sf::st_transform(region, sf::st_crs(geography)))
}
selected <- inner_join(geography, green, by = "GEOID")
p3 <- ggplot() +
  geom_sf(data = sf::st_union(geography), fill = "#f0f0f0", color = "#777777", linewidth = 0.3) +
  geom_sf(data = selected, aes(fill = pct_green_of_developed * 100), color = "#303030", linewidth = 0.2) +
  scale_fill_gradientn(colors = c("#d95f02", "#f7e225", "#1a9850"), limits = c(0, 100)) +
  labs(title = "Greenness extraction footprint", subtitle = paste(nrow(selected), "sampled block groups | scope", scope),
    fill = "Green (%)", caption = "Only highlighted block groups were processed. Unshaded land has not been measured. Map is not filtered by supplier GPCD.")
ggsave(file.path(output_dir, "sample_map.png"), p3, width = 8, height = 8, dpi = 180)
  cat("Suppliers with usable records:", nrow(df_all), "\nSuppliers above GPCD maximum:", excluded_high_gpcd,
      "\nSuppliers plotted:", nrow(df), "\nSaved three figures to", output_dir, "\n")
