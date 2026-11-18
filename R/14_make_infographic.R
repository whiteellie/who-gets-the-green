source(here::here("R", "00_config.R"))
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

ink <- "#1c2b25"
paper <- "#f7f8f3"
green_dark <- "#1f6548"
panel_theme <- theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    plot.title = element_text(color = ink, face = "bold", size = 12),
    plot.subtitle = element_text(color = "#53645b", size = 9),
    plot.margin = margin(9, 12, 8, 12),
    panel.grid = element_blank()
  )

greenness_files <- list.files(
  dir_raw,
  pattern = "^greenness_bg_hr_[a-z0-9_]+[.]csv$",
  full.names = TRUE
)
greenness_files <- greenness_files[!grepl("_errors[.]csv$", greenness_files)]
regions <- read_sf(file.path(dir_raw, "hydrologic_regions.gpkg"), quiet = TRUE) |>
  st_transform(3310) |>
  st_make_valid()
region_slugs <- gsub("[^a-z0-9]+", "_", tolower(regions$HR_NAME))
block_groups <- read_sf(file.path(dir_raw, "acs_bg.gpkg"), quiet = TRUE) |>
  st_transform(3310) |>
  st_make_valid()

map_parts <- lapply(greenness_files, function(path) {
  scope <- sub("^greenness_bg_(hr_[a-z0-9_]+)[.]csv$", "\\1", basename(path))
  region_slug <- sub("^hr_", "", scope)
  region_index <- which(region_slugs == region_slug)
  stopifnot(length(region_index) == 1L)
  green <- read_csv(
    path,
    col_select = c(GEOID, pct_green_of_developed),
    col_types = cols(GEOID = "c", pct_green_of_developed = col_double()),
    show_col_types = FALSE
  )
  stopifnot(!anyDuplicated(green$GEOID))
  region_bg <- block_groups |>
    semi_join(green, by = "GEOID") |>
    st_intersection(regions[region_index, "HR_NAME"]) |>
    left_join(green, by = "GEOID")
  stopifnot(nrow(region_bg) == nrow(green))
  region_bg
})
greenness_map <- do.call(rbind, map_parts)
regions_map <- regions
stopifnot(
  nrow(greenness_map) == sum(vapply(map_parts, nrow, integer(1))),
  all(is.na(greenness_map$pct_green_of_developed) |
    dplyr::between(greenness_map$pct_green_of_developed, 0, 1))
)

header <- ggplot() +
  annotate("text", x = 0, y = 0.82, label = "WHO GETS THE GREEN?", hjust = 0,
    color = ink, fontface = "bold", size = 8.5) +
  annotate("text", x = 0, y = 0.22,
    label = "A California snapshot of developed-land greenness and what helps explain it",
    hjust = 0, color = "#53645b", size = 4) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = paper, color = NA), plot.margin = margin(8, 14, 2, 14))

map_panel <- ggplot() +
  geom_sf(data = greenness_map, aes(fill = pct_green_of_developed), color = NA) +
  geom_sf(data = regions_map, fill = NA, color = "#53645b", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = c("#f2e7ad", "#d8ed9e", "#78c679", "#238443", "#005a32"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    na.value = "#cbd0ca",
    name = "Green share"
  ) +
  labs(
    title = "1. Where is it green?",
    subtitle = "Developed-land green share; gray = no data"
  ) +
  coord_sf(datum = NA, expand = FALSE) +
  panel_theme +
  theme(
    axis.text = element_blank(), axis.title = element_blank(),
    legend.position = "bottom", legend.key.width = grid::unit(1.3, "cm"),
    plot.background = element_rect(fill = paper, color = NA)
  )

powerpoint_panel_dir <- file.path(dir_figures, "powerpoint_panels")
dir.create(powerpoint_panel_dir, recursive = TRUE, showWarnings = FALSE)
greenness_powerpoint_map <- map_panel +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  guides(fill = "none") +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 0, 6, 6)
  )
green_legend_data <- data.frame(x = 0.1, green_share = seq(0.005, 0.995, length.out = 200))
green_legend_labels <- data.frame(
  x = 0.23,
  green_share = c(0.94, 0.73, 0.5, 0.27, 0.06),
  label = scales::label_percent(accuracy = 1)(seq(1, 0, by = -0.25))
)
green_legend <- ggplot(green_legend_data, aes(x, green_share, fill = green_share)) +
  geom_tile(width = 0.16, height = 0.006) +
  geom_text(
    data = green_legend_labels,
    aes(x, green_share, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 18 / ggplot2::.pt,
    color = ink
  ) +
  scale_fill_gradientn(
    colors = c("#f2e7ad", "#d8ed9e", "#78c679", "#238443", "#005a32"),
    limits = c(0, 1), guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 0.9), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1), breaks = NULL, expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 0)
  )
greenness_powerpoint_plot <- (greenness_powerpoint_map | green_legend) +
  plot_layout(widths = c(0.76, 0.24))
ggsave(
  file.path(powerpoint_panel_dir, "01_california_greenness_map.png"),
  greenness_powerpoint_plot,
  width = 6.5,
  height = 6.5,
  dpi = 300,
  bg = paper
)

response_panel <- ggplot() +
  annotate("text", x = 0.04, y = 0.96, label = "THE RESPONSE",
    hjust = 0, vjust = 1, color = ink, fontface = "bold", size = 4.1) +
  annotate("text", x = 0.04, y = 0.79,
    label = "Fraction of valid developed\npixels classified green",
    hjust = 0, vjust = 1, color = green_dark, fontface = "bold", size = 4, lineheight = 1.05) +
  annotate("text", x = 0.04, y = 0.55,
    label = "NDVI = (NIR − red) / (NIR + red)\nGreen pixel: NDVI > 0.30",
    hjust = 0, vjust = 1, color = ink, size = 2.8, lineheight = 1.12) +
  annotate("text", x = 0.04, y = 0.32,
    label = "NLCD 2021 developed mask\nSentinel-2 summer 2022",
    hjust = 0, vjust = 1, color = ink, size = 2.8, lineheight = 1.12) +
  annotate("text", x = 0.04, y = 0.08,
    label = "PWS outcome: population-weighted\nmean over sampled overlaps.",
    hjust = 0, vjust = 0, color = "#53645b", size = 2.7, lineheight = 1.08) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = paper, color = NA), plot.margin = margin(12, 14, 12, 10))

predictor_steps <- data.frame(
  x = seq_len(5),
  model = paste0("M", seq_len(5)),
  predictor = c(
    "Outdoor use + region\nHow much outdoor water\nper person?",
    "Climate\nHow hot and dry?\nETo + precipitation",
    "Development\nHow densely built is\nthe community?",
    "Landscape space\nHow much developed land\nis non-impervious?",
    "Income proxy\nDoes median income help\nexplain what's left?"
  )
)
predictor_panel <- ggplot() +
  geom_segment(
    aes(x = 1.12, xend = 4.82, y = 0.55, yend = 0.55),
    color = "#9ba99f", linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.13, "in"), type = "closed")
  ) +
  geom_point(data = predictor_steps, aes(x, 0.55), color = green_dark, size = 5) +
  geom_text(data = predictor_steps, aes(x, 0.82, label = model),
    color = "#53645b", fontface = "bold", size = 3.1) +
  geom_text(data = predictor_steps, aes(x, 0.23, label = predictor),
    color = ink, size = 3, lineheight = 1.1) +
  scale_x_continuous(breaks = NULL, limits = c(0.55, 5.45)) +
  scale_y_continuous(limits = c(0, 1), breaks = NULL) +
  labs(
    title = "2. What explains the greenness?",
    subtitle = "Terms accumulate; pooled fits include region effects."
  ) +
  panel_theme +
  theme(
    axis.title = element_blank(),
    plot.background = element_rect(fill = paper, color = NA)
  )

incremental <- read_csv(
  file.path(here::here("outputs", "models"), "incremental_cv_r_squared.csv"),
  show_col_types = FALSE
) |>
  filter(grepl("^M[1-5]:", as.character(model))) |>
  arrange(match(as.character(model), c(
    "M1: use + region, <=132", "M2: add climate, <=132", "M3: add development, <=132",
    "M4: add landscape space, <=132", "M5: add income, <=132"
  ))) |>
  mutate(
    step = factor(paste0("M", seq_len(n())), levels = paste0("M", 1:5)),
    model_label = factor(paste0("M", seq_len(n()), "\n", c("Use", "Climate", "Developed", "Space", "Income")),
      levels = paste0("M", 1:5, "\n", c("Use", "Climate", "Developed", "Space", "Income"))),
    bar_label = if_else(is.na(incremental_r_squared),
      sprintf("%.3f\nbaseline", r_squared),
      sprintf("%.3f\n%+.3f", r_squared, incremental_r_squared))
  )
stopifnot(nrow(incremental) == 5L, all(incremental$observations == incremental$observations[[1]]))

residuals <- read_csv(
  file.path(here::here("outputs", "models"), "all_regions_m5_cv_residuals.csv"),
  show_col_types = FALSE
)
spatial_map_path <- file.path(here::here("outputs", "models"), "all_regions_m5_spatial_residuals.gpkg")
if (!file.exists(spatial_map_path)) stop("Run R/15_spatial_residual_analysis.R before building the infographic.")
spatial_residuals <- read_sf(spatial_map_path, quiet = TRUE)
anomaly_labels <- bind_rows(
  filter(residuals, residual_pct > 0) |>
    arrange(desc(residual_pct)) |>
    slice_head(n = 3) |>
    mutate(point_label = paste0("G", row_number()), label_y = c(0.74, 0.64, 0.54)),
  filter(residuals, residual_pct < 0) |>
    arrange(residual_pct) |>
    slice_head(n = 3) |>
    mutate(point_label = paste0("L", row_number()), label_y = c(0.31, 0.21, 0.11))
)
stopifnot(nrow(anomaly_labels) == 6L, all(!is.na(anomaly_labels$supplier_name)))

r_squared_panel <- ggplot(incremental, aes(model_label, r_squared, fill = step)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = bar_label), vjust = -0.3, color = ink, size = 3.05, lineheight = 0.95) +
  scale_fill_manual(values = c("#7f9587", "#668c70", "#69965f", "#3d8c62", "#1f6548"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2),
    labels = scales::label_number(accuracy = 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "3. How much can we explain?",
    subtitle = "Cumulative held-out R²; change from prior step",
    x = NULL,
    y = "Five-fold CV R²"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(color = ink, size = 8.5),
    axis.title.y = element_text(color = "#53645b", size = 8),
    plot.background = element_rect(fill = paper, color = NA)
  )

residual_plot <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#7b8980", linewidth = 0.45) +
  geom_point(
    data = filter(residuals, !anomaly_flag),
    aes(observed_pct, predicted_pct), color = "#bac2bc", alpha = 0.45, size = 0.9
  ) +
  geom_point(
    data = filter(residuals, anomaly_flag),
    aes(observed_pct, predicted_pct, fill = residual_class),
    shape = 21, color = ink, stroke = 0.25, size = 2.2
  ) +
  geom_text(
    data = anomaly_labels,
    aes(observed_pct, predicted_pct, label = point_label),
    color = "white", fontface = "bold", size = 1.8
  ) +
  scale_fill_manual(
    values = c("Less green than expected" = "#c45b42", "Greener than expected" = "#258778"),
    breaks = c("Greener than expected", "Less green than expected"),
    labels = c("Greener", "Less green"),
    name = "Large residual"
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = 0.02)) +
  coord_equal() +
  labs(
    title = "4. How did we do?",
    subtitle = paste0(
      "Final M5 out-of-fold • ", sum(residuals$anomaly_flag), " large residuals (|e| ≥ ",
      sprintf("%.1f", unique(residuals$anomaly_cutoff_pct)), " pp)"
    ),
    x = "Observed green (%)",
    y = "Predicted green (%)",
    caption = "Point IDs map to panel 06; cutoff is a review screen, not an outlier test."
  ) +
  panel_theme +
  theme(
    axis.text = element_text(size = 7, color = ink),
    axis.title = element_text(size = 8, color = "#53645b"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 8),
    plot.caption = element_text(size = 7, hjust = 0),
    plot.background = element_rect(fill = paper, color = NA)
  )

residual_panel <- residual_plot

spatial_residual_map_panel <- ggplot() +
  geom_sf(data = spatial_residuals, aes(fill = residual_band), color = NA) +
  geom_sf(data = regions, fill = NA, color = "#ffffff", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "≤ −20 pp" = "#bd4f39",
      "−20 to −10 pp" = "#e5a08c",
      "−10 to +10 pp" = "#e7ebdf",
      "+10 to +20 pp" = "#84c4ad",
      "≥ +20 pp" = "#167565"
    ),
    drop = FALSE,
    name = "Observed − predicted"
  ) +
  guides(fill = guide_legend(ncol = 1, keyheight = grid::unit(0.28, "cm"))) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "5. Where are we off?",
    subtitle = "M5 residual: observed minus predicted greenness"
  ) +
  panel_theme +
  theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank(), legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 5.5),
    plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8),
    plot.background = element_rect(fill = paper, color = NA)
  )

residual_band_key <- data.frame(
  y = seq(5, 1),
  band = factor(c("≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp"),
    levels = c("≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp")),
  label = c("≤ −20", "−20 to −10", "−10 to +10", "+10 to +20", "≥ +20")
)
residual_band_legend <- ggplot(residual_band_key, aes(0.16, y, fill = band)) +
  geom_tile(width = 0.24, height = 0.82) +
  geom_text(aes(x = 0.34, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 2.25, color = ink) +
  scale_fill_manual(
    values = c(
      "≤ −20 pp" = "#bd4f39", "−20 to −10 pp" = "#e5a08c", "−10 to +10 pp" = "#e7ebdf",
      "+10 to +20 pp" = "#84c4ad", "≥ +20 pp" = "#167565"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1.3), breaks = NULL) +
  scale_y_continuous(limits = c(0.5, 5.5), breaks = NULL) +
  theme_void() +
  theme(plot.background = element_rect(fill = paper, color = NA), plot.margin = margin(8, 8, 8, 0))
spatial_residual_panel <- spatial_residual_map_panel | residual_band_legend +
  plot_layout(widths = c(0.78, 0.22))

standout_labels <- bind_rows(
  filter(anomaly_labels, residual_pct > 0) |>
    mutate(label_y = c(0.76, 0.64, 0.52)),
  filter(anomaly_labels, residual_pct < 0) |>
    mutate(label_y = c(0.30, 0.18, 0.06))
)
standouts_panel <- ggplot() +
  annotate("text", x = 0.04, y = 0.96, label = "6. Which communities defy expectations?",
    hjust = 0, vjust = 1, color = ink, fontface = "bold", size = 4) +
  annotate("text", x = 0.04, y = 0.85, label = "GREENER THAN EXPECTED",
    hjust = 0, color = "#17685f", fontface = "bold", size = 3) +
  geom_text(
    data = filter(standout_labels, residual_pct > 0),
    aes(x = 0.04, y = label_y, label = paste0(point_label, "  ", supplier_name)),
    hjust = 0, color = ink, size = 2.7
  ) +
  geom_text(
    data = filter(standout_labels, residual_pct > 0),
    aes(x = 0.96, y = label_y, label = sprintf("%+.1f pp", residual_pct)),
    hjust = 1, color = "#17685f", fontface = "bold", size = 2.6
  ) +
  annotate("text", x = 0.04, y = 0.43, label = "LESS GREEN THAN EXPECTED",
    hjust = 0, color = "#a84331", fontface = "bold", size = 3) +
  geom_text(
    data = filter(standout_labels, residual_pct < 0),
    aes(x = 0.04, y = label_y, label = paste0(point_label, "  ", supplier_name)),
    hjust = 0, color = ink, size = 2.7
  ) +
  geom_text(
    data = filter(standout_labels, residual_pct < 0),
    aes(x = 0.96, y = label_y, label = sprintf("%+.1f pp", residual_pct)),
    hjust = 1, color = "#a84331", fontface = "bold", size = 2.6
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = "#eef1e8", color = NA), plot.margin = margin(15, 14, 8, 14))

footer <- ggplot() +
  annotate("text", x = 0, y = 0.5,
    label = paste0("1,282 supplier-region records  •  1,234 PWSs  •  M5 CV R² = ",
      sprintf("%.3f", incremental$r_squared[[5]]),
      "  •  Exploratory associations, not causal estimates or efficiency rankings"),
    hjust = 0, color = "#53645b", size = 3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = paper, color = NA), plot.margin = margin(3, 14, 6, 14))

panel_frame <- theme(
  plot.background = element_rect(fill = paper, color = "#d6e0d9", linewidth = 0.5),
  plot.margin = margin(10, 12, 10, 12)
)
top_panels <- (
  (patchwork::wrap_elements(full = (map_panel | response_panel)) | predictor_panel | r_squared_panel) +
    plot_layout(widths = c(1, 1, 1))
) & panel_frame
bottom_panels <- (
  (residual_panel | patchwork::wrap_elements(full = spatial_residual_panel) | standouts_panel) +
    plot_layout(widths = c(1, 1, 1))
) & panel_frame

infographic <- (
  header /
    top_panels /
    bottom_panels /
    footer
) +
  plot_layout(heights = c(0.17, 1, 1, 0.08)) &
  theme(plot.background = element_rect(fill = paper, color = NA))

output_path <- file.path(dir_figures, "california_greenness_infographic.png")
ggsave(output_path, infographic, width = 18, height = 13.5, dpi = 220, bg = paper)
cat(
  "Mapped block-group pieces:", nrow(greenness_map),
  "\nNo-valid-green values:", sum(!is.finite(greenness_map$pct_green_of_developed)),
  "\nSaved infographic to", output_path, "\n"
)