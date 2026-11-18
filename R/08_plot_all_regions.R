source(here::here("R", "00_config.R"))
library(dplyr)
library(readr)
library(ggplot2)

max_outdoor_gpcd <- 132
region_scopes <- list.dirs(dir_processed, recursive = FALSE, full.names = FALSE) |>
  grep(pattern = "^hr_", value = TRUE) |>
  sort()
stopifnot(length(region_scopes) == 10L)

region_data <- lapply(region_scopes, function(scope) {
  region_name <- tools::toTitleCase(gsub("_", " ", sub("^hr_", "", scope)))
  read_csv(
    file.path(dir_processed, scope, "pws_green_water_income.csv"),
    show_col_types = FALSE
  ) |>
    mutate(region = region_name)
})

df_all <- bind_rows(region_data) |>
  filter(
    water_data_complete,
    is.finite(outdoor_gpcd),
    is.finite(pct_green_sampled),
    is.finite(income_proxy)
  )
excluded_high_gpcd <- sum(df_all$outdoor_gpcd > max_outdoor_gpcd)
df <- df_all |>
  filter(outdoor_gpcd <= max_outdoor_gpcd) |>
  mutate(region = factor(region, levels = unique(region)))
stopifnot(nrow(df) > 0L, n_distinct(df$region) == 10L)

caption <- paste(
  "Each point is a supplier-region record; suppliers spanning regions may appear more than once.",
  "Water use is supplier-wide and is not apportioned to regions.",
  paste0(excluded_high_gpcd, " records above ", max_outdoor_gpcd, " GPCD excluded; source tables unchanged."),
  "Five Central Coast, two South Lahontan, and five Colorado River block groups remain unprocessed.",
  "Water and imagery: 2022. ACS: 2019-2023. No climate adjustment; not an efficiency ranking.",
  sep = "\n"
)

figure_dir <- file.path(dir_figures, "all_regions_gpcd_max_132")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
region_colors <- setNames(scales::hue_pal()(nlevels(df$region)), levels(df$region))

plot <- ggplot(df, aes(outdoor_gpcd, pct_green_sampled, fill = region, size = income_proxy)) +
  geom_point(shape = 21, color = "black", stroke = 0.45, alpha = 0.72) +
  scale_fill_manual(values = region_colors, name = "Hydrologic region") +
  scale_size_continuous(name = "Income proxy", labels = scales::label_dollar(), range = c(1.5, 8)) +
  scale_x_continuous(limits = c(0, max_outdoor_gpcd), breaks = seq(0, max_outdoor_gpcd, 20)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Estimated outdoor residential use (gallons/person/day)",
    y = "Sampled developed-land greenness (%)",
    title = "Greenness and outdoor use across hydrologic regions",
    subtitle = paste(nrow(df), "supplier-region records across", n_distinct(df$region), "regions"),
    caption = caption
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(
  file.path(figure_dir, "all_regions_green_vs_outdoor_use.png"),
  plot,
  width = 11,
  height = 8,
  dpi = 180
)

cat(
  "Supplier-region records with usable data:", nrow(df_all),
  "\nRecords above GPCD maximum:", excluded_high_gpcd,
  "\nSupplier-region records plotted:", nrow(df),
  "\nSaved combined figure to", figure_dir, "\n"
)