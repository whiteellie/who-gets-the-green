source(here::here("R", "00_config.R"))
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

set.seed(2022)
permutations <- 9999L
neighbors_k <- 8L

residuals <- read_csv(
  file.path(here::here("outputs", "models"), "all_regions_m5_cv_residuals.csv"),
  show_col_types = FALSE
)
stopifnot(!anyDuplicated(residuals[c("pwsid", "region")]))

pws <- read_sf(file.path(dir_raw, "pws_boundaries.gpkg"), quiet = TRUE) |>
  filter(BOUNDARY_TYPE == "Water Service Area", !is.na(SABL_PWSID)) |>
  st_transform(3310) |>
  st_make_valid() |>
  transmute(pwsid = SABL_PWSID) |>
  group_by(pwsid) |>
  summarise(do_union = TRUE, .groups = "drop")
regions <- read_sf(file.path(dir_raw, "hydrologic_regions.gpkg"), quiet = TRUE) |>
  st_transform(3310) |>
  st_make_valid()
state_outline <- st_sf(geometry = st_union(st_geometry(regions)))
stopifnot(all(residuals$pwsid %in% pws$pwsid))

spatial_parts <- lapply(sort(unique(residuals$region)), function(region_name) {
  region <- filter(regions, HR_NAME == region_name)
  records <- filter(residuals, region == region_name)
  clipped <- pws |>
    semi_join(records, by = "pwsid") |>
    st_intersection(region |> select(HR_NAME)) |>
    left_join(records, by = "pwsid")
  stopifnot(
    nrow(clipped) == nrow(records),
    !anyDuplicated(clipped[c("pwsid", "region")]),
    all(is.finite(clipped$residual_pct)),
    all(st_is_valid(clipped))
  )
  clipped
})
spatial_residuals <- do.call(rbind, spatial_parts)
stopifnot(nrow(spatial_residuals) == nrow(residuals))

neighbor_graph <- function(data, k = 8L) {
  count <- nrow(data)
  if (count < 3L) return(NULL)
  point_coordinates <- st_coordinates(st_point_on_surface(data))[, 1:2, drop = FALSE]
  distances <- as.matrix(stats::dist(point_coordinates))
  diag(distances) <- Inf
  k <- min(k, count - 1L)
  neighbors <- matrix(FALSE, nrow = count, ncol = count)
  for (row in seq_len(count)) {
    neighbors[row, order(distances[row, ])[seq_len(k)]] <- TRUE
  }
  neighbors <- neighbors | t(neighbors)
  diag(neighbors) <- FALSE
  edge <- which(neighbors, arr.ind = TRUE)
  degree <- rowSums(neighbors)
  list(
    n = count,
    from = edge[, 1],
    to = edge[, 2],
    weight = 1 / degree[edge[, 1]],
    residual = data$residual_pct - mean(data$residual_pct),
    k = k
  )
}

moran_statistic <- function(values, graph) {
  n <- length(values)
  centered <- values - mean(values)
  n / sum(graph$weight) *
    sum(centered[graph$from] * centered[graph$to] * graph$weight) /
    sum(centered^2)
}

region_graphs <- lapply(sort(unique(spatial_residuals$region)), function(region_name) {
  region_data <- filter(spatial_residuals, region == region_name)
  graph <- neighbor_graph(region_data, neighbors_k)
  list(region = region_name, graph = graph)
})
region_graphs <- Filter(function(item) !is.null(item$graph), region_graphs)
stopifnot(length(region_graphs) == n_distinct(spatial_residuals$region))

region_tests <- lapply(region_graphs, function(item) {
  graph <- item$graph
  observed <- moran_statistic(graph$residual, graph)
  simulated <- replicate(permutations, moran_statistic(sample(graph$residual), graph))
  center <- mean(simulated)
  data.frame(
    region = item$region,
    observations = graph$n,
    k_neighbors = graph$k,
    morans_i = observed,
    permutation_p = (1 + sum(abs(simulated - center) >= abs(observed - center))) / (permutations + 1),
    permutations = permutations
  )
}) |>
  bind_rows() |>
  mutate(permutation_p_bh = p.adjust(permutation_p, method = "BH")) |>
  arrange(permutation_p)

offsets <- cumsum(c(0L, vapply(region_graphs, function(item) item$graph$n, integer(1))))
pooled_graph <- list(
  n = sum(vapply(region_graphs, function(item) item$graph$n, integer(1))),
  from = unlist(Map(function(item, offset) item$graph$from + offset, region_graphs, offsets[-length(offsets)])),
  to = unlist(Map(function(item, offset) item$graph$to + offset, region_graphs, offsets[-length(offsets)])),
  weight = unlist(lapply(region_graphs, function(item) item$graph$weight))
)
pooled_residual <- unlist(lapply(region_graphs, function(item) item$graph$residual))
pooled_graph$residual <- pooled_residual
pooled_observed <- moran_statistic(pooled_residual, pooled_graph)
pooled_simulated <- replicate(permutations, {
  permuted <- unlist(lapply(region_graphs, function(item) sample(item$graph$residual)))
  moran_statistic(permuted, pooled_graph)
})
pooled_center <- mean(pooled_simulated)
pooled_test <- data.frame(
  region = "California, stratified within region",
  observations = pooled_graph$n,
  k_neighbors = neighbors_k,
  morans_i = pooled_observed,
  permutation_p = (1 + sum(abs(pooled_simulated - pooled_center) >= abs(pooled_observed - pooled_center))) /
    (permutations + 1),
  permutations = permutations,
  permutation_p_bh = NA_real_
)
moran_results <- bind_rows(pooled_test, region_tests)
spatial_residuals <- mutate(
  spatial_residuals,
  residual_band = cut(
    residual_pct,
    breaks = c(-Inf, -20, -10, 10, 20, Inf),
    labels = c("≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp"),
    right = TRUE
  ),
  residual_band = factor(residual_band, levels = c(
    "≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp"
  ))
)

model_dir <- here::here("outputs", "models")
figure_dir <- file.path(dir_figures, "all_regions_model")
powerpoint_panel_dir <- file.path(dir_figures, "powerpoint_panels")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(powerpoint_panel_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(moran_results, file.path(model_dir, "all_regions_m5_residual_morans_i.csv"))
st_write(
  spatial_residuals,
  file.path(model_dir, "all_regions_m5_spatial_residuals.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

ca_map_standouts <- bind_rows(
  filter(spatial_residuals, residual_pct > 0) |>
    arrange(desc(residual_pct)) |>
    slice_head(n = 3) |>
    mutate(map_id = paste0("G", row_number())),
  filter(spatial_residuals, residual_pct < 0) |>
    arrange(residual_pct) |>
    slice_head(n = 3) |>
    mutate(map_id = paste0("L", row_number()))
)
ca_map_standout_points <- st_point_on_surface(ca_map_standouts)
ca_map_coordinates <- st_coordinates(ca_map_standout_points)
ca_map_standout_labels <- st_drop_geometry(ca_map_standout_points) |>
  mutate(label_x = ca_map_coordinates[, 1], label_y = ca_map_coordinates[, 2])
stopifnot(nrow(ca_map_standout_labels) == 6L)

residual_map <- ggplot() +
  geom_sf(data = state_outline, fill = "#e3e6e3", color = "#7b8980", linewidth = 0.25) +
  geom_sf(data = spatial_residuals, aes(fill = residual_band), color = NA) +
  geom_sf(data = regions, fill = NA, color = "#ffffff", linewidth = 0.22) +
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
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "Greenness residuals cluster geographically",
    subtitle = paste(nrow(spatial_residuals), "PWS-region records • M5 PWS-grouped out-of-fold residuals"),
    caption = paste(
      "Service areas are clipped to HR boundaries; color represents the sampled-overlap residual, not every parcel in the area.",
      "Green = observed greener than predicted; rust = observed less green. Colors group residuals into percentage-point bands.",
      sep = "\n"
    ),
    fill = "Residual (pp)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  file.path(figure_dir, "m5_residual_map.png"),
  residual_map,
  width = 8,
  height = 11,
  dpi = 180
)
residual_powerpoint_plot <- residual_map +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  guides(fill = "none") +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )
residual_legend_data <- data.frame(
  y = 5:1,
  residual_band = factor(
    c("≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp"),
    levels = c("≤ −20 pp", "−20 to −10 pp", "−10 to +10 pp", "+10 to +20 pp", "≥ +20 pp")
  ),
  label = c("≤ −20", "−20 to −10", "−10 to +10", "+10 to +20", "≥ +20")
)
residual_legend <- ggplot(residual_legend_data, aes(0.12, y, fill = residual_band)) +
  geom_tile(width = 0.16, height = 0.82) +
  geom_text(aes(x = 0.36, y = y, label = label), inherit.aes = FALSE,
    hjust = 0, size = 18 / ggplot2::.pt, color = "#1c2b25") +
  scale_fill_manual(
    values = c(
      "≤ −20 pp" = "#bd4f39", "−20 to −10 pp" = "#e5a08c", "−10 to +10 pp" = "#e7ebdf",
      "+10 to +20 pp" = "#84c4ad", "≥ +20 pp" = "#167565"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1.7), breaks = NULL) +
  scale_y_continuous(limits = c(0.5, 5.5), breaks = NULL) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 0)
  )
residual_powerpoint_panel <- (residual_powerpoint_plot | residual_legend) +
  patchwork::plot_layout(widths = c(0.74, 0.26))
ggsave(
  file.path(powerpoint_panel_dir, "02_final_model_residual_map.png"),
  residual_powerpoint_panel,
  width = 6.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

standout_region_map <- ggplot() +
  geom_sf(data = regions, fill = "#e3e8e3", color = "#ffffff", linewidth = 0.65) +
  geom_sf(data = state_outline, fill = NA, color = "#53645b", linewidth = 0.7) +
  geom_point(
    data = ca_map_standout_labels,
    aes(label_x, label_y, color = residual_class),
    size = 3.2
  ) +
  ggrepel::geom_label_repel(
    data = ca_map_standout_labels,
    aes(label_x, label_y, label = map_id, color = residual_class),
    seed = 2022,
    size = 18 / ggplot2::.pt,
    box.padding = 0.55,
    point.padding = 0.5,
    label.padding = grid::unit(0.22, "lines"),
    label.r = grid::unit(0.08, "lines"),
    label.size = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 10,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Less green than expected" = "#a84331", "Greener than expected" = "#17685f"),
    guide = "none"
  ) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave(
  file.path(powerpoint_panel_dir, "05_standout_hydrologic_regions.png"),
  standout_region_map,
  width = 6.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

standout_systems <- bind_rows(
  filter(spatial_residuals, residual_pct > 0) |>
    arrange(desc(residual_pct)) |>
    slice_head(n = 3) |>
    mutate(direction = "Greener than expected"),
  filter(spatial_residuals, residual_pct < 0) |>
    arrange(residual_pct) |>
    slice_head(n = 3) |>
    mutate(direction = "Less green than expected")
)
stopifnot(nrow(standout_systems) == 6L, all(!is.na(standout_systems$supplier_name)))

plot_standout_map <- function(direction, output_name, highlight, label_color, id_prefix) {
  selected <- filter(standout_systems, direction == .env$direction) |>
    arrange(desc(abs(residual_pct))) |>
    mutate(map_id = paste0(id_prefix, row_number()))
  selected_regions <- unique(selected$region)

  local_maps <- lapply(selected_regions, function(region_name) {
    region_selected <- filter(selected, region == region_name)
    context <- filter(spatial_residuals, region == region_name)
    context_region <- filter(regions, HR_NAME == region_name)
    points <- st_point_on_surface(region_selected)
    coordinates <- st_coordinates(points)
    point_labels <- st_drop_geometry(points) |>
      mutate(label_x = coordinates[, 1], label_y = coordinates[, 2])
    bounds <- st_bbox(region_selected)
    span <- max(bounds["xmax"] - bounds["xmin"], bounds["ymax"] - bounds["ymin"], 20000) * 2.2
    center_x <- (bounds["xmin"] + bounds["xmax"]) / 2
    center_y <- (bounds["ymin"] + bounds["ymax"]) / 2
    x_limits <- c(center_x - span / 2, center_x + span / 2)
    y_limits <- c(center_y - span / 2, center_y + span / 2)

    ggplot() +
      geom_sf(data = context_region, fill = "#eef0ed", color = "#79867d", linewidth = 0.35) +
      geom_sf(data = context, fill = "#e2e7e2", color = "#b8c0b9", linewidth = 0.12) +
      geom_sf(data = region_selected, fill = highlight, color = label_color, linewidth = 0.55) +
      ggrepel::geom_label_repel(
        data = point_labels,
        aes(label_x, label_y, label = map_id),
        seed = 2022,
        size = 14 / ggplot2::.pt,
        box.padding = 0.35,
        point.padding = 0.2,
        min.segment.length = 0,
        max.overlaps = Inf,
        max.time = 5,
        label.padding = grid::unit(0.14, "lines"),
        label.r = grid::unit(0.08, "lines"),
        label.size = 0.2,
        fill = "white",
        color = label_color,
        fontface = "bold"
      ) +
      coord_sf(xlim = x_limits, ylim = y_limits, datum = NA, expand = FALSE) +
      labs(title = region_name) +
      theme_void(base_size = 14, base_family = "sans") +
      theme(
        plot.title = element_text(face = "bold", color = "#1c2b25", size = 16),
        plot.title.position = "plot",
        plot.margin = margin(8, 8, 8, 8),
        plot.background = element_rect(fill = "white", color = "#d6e0d9", linewidth = 0.4)
      )
  })

  key_data <- st_drop_geometry(selected) |>
    transmute(map_id, supplier_name, residual_pct) |>
    mutate(
      y = rev(seq_len(n())),
      key_text = paste0(map_id, "  ", supplier_name, "\n     ", sprintf("%+.1f pp", residual_pct))
    )
  key_panel <- ggplot(key_data) +
    geom_text(aes(x = 0, y, label = key_text), hjust = 0, size = 12 / ggplot2::.pt,
      color = label_color, lineheight = 1.05) +
    scale_x_continuous(limits = c(0, 1), breaks = NULL) +
    scale_y_continuous(limits = c(0.5, nrow(key_data) + 0.5), breaks = NULL) +
    labs(title = "Highlighted systems") +
    theme_void(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", color = "#1c2b25", size = 15),
      plot.title.position = "plot",
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "#f4f6f2", color = "#d6e0d9", linewidth = 0.4)
    )
  map_panels <- c(local_maps, list(key_panel))
  map_widths <- c(rep(1, length(local_maps)), 1.35)
  figure <- patchwork::wrap_plots(map_panels, ncol = length(map_panels), widths = map_widths) +
    patchwork::plot_annotation(
      title = direction,
      theme = theme(plot.title = element_text(face = "bold", size = 22, color = "#1c2b25"))
    )
  ggsave(file.path(powerpoint_panel_dir, output_name), figure, width = 13, height = 7.5, dpi = 300, bg = "white")
}

plot_standout_map("Greener than expected", "05_top_greener_than_expected_map.png", "#258778", "#145e55", "G")
plot_standout_map("Less green than expected", "06_top_less_green_than_expected_map.png", "#c45b42", "#903b2b", "L")

cat("Mapped PWS-region records:", nrow(spatial_residuals), "\n")
cat("Within-region 8-nearest-neighbor Moran tests (BH-adjusted):\n")
print(moran_results, row.names = FALSE)
cat("Saved map to", file.path(figure_dir, "m5_residual_map.png"), "\n")