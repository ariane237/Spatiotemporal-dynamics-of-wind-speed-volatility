# ============================================================
# Script: 01_load_prepare_and_preprocess_wind_data.R
# Purpose:
#   Load Agrimonia wind-speed data, create descriptive summaries,
#   generate illustrative station plots, and construct STL + AR(1)
#   preprocessed panels for ws10 and ws100.
#
# Author: Ariane Meli Chrisko
# ============================================================

# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(data.table)
library(ggplot2)
library(moments)
library(zoo)
library(pheatmap)
library(viridis)
library(sf)
library(dplyr)
library(lubridate)
library(leaflet)
library(htmltools)
library(htmlwidgets)
library(webshot2)
library(magick)
library(patchwork)

# ------------------------------------------------------------
# 2. Define file paths
# ------------------------------------------------------------

agrimonia_file <- "data/raw/Agrimonia_Dataset_v_3_0_0.csv"
nuts_dir <- "data/raw/NUTS"
output_dir <- "results"

dir.create(nuts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Load and prepare Agrimonia wind data
# ------------------------------------------------------------

dt <- fread(agrimonia_file)

dt[, Date := as.Date(Time, format = "%d/%m/%Y")]

wind <- dt[, .(
  IDStations,
  Latitude,
  Longitude,
  Date,
  ws10_mean  = WE_wind_speed_10m_mean,
  ws10_max   = WE_wind_speed_10m_max,
  ws100_mean = WE_wind_speed_100m_mean,
  ws100_max  = WE_wind_speed_100m_max
)]

wind[, IDStations := as.character(IDStations)]

# ------------------------------------------------------------
# 4. Load NUTS shapefile and extract Lombardy
# ------------------------------------------------------------

nuts_file <- file.path(nuts_dir, "NUTS_RG_01M_2021_4326.shp")
nuts2 <- st_read(nuts_file, quiet = TRUE)

lombardy <- nuts2[nuts2$NUTS_ID == "ITC4", ]
lombardy_sf <- st_transform(lombardy, 4326)

# ------------------------------------------------------------
# 5. Create leaflet map of monitoring stations
# ------------------------------------------------------------

stations_unique <- unique(wind[, .(IDStations, Longitude, Latitude)])
stations_sf <- st_as_sf(
  stations_unique,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

bb <- st_bbox(lombardy_sf)
pad <- 0.03

xmin <- as.numeric(bb["xmin"] + pad * (bb["xmax"] - bb["xmin"]))
xmax <- as.numeric(bb["xmax"] - pad * (bb["xmax"] - bb["xmin"]))
ymin <- as.numeric(bb["ymin"] + pad * (bb["ymax"] - bb["ymin"]))
ymax <- as.numeric(bb["ymax"] - pad * (bb["ymax"] - bb["ymin"]))

map_lombardy <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
  addProviderTiles(providers$CartoDB.Voyager) %>%
  fitBounds(xmin, ymin, xmax, ymax) %>%
  addPolygons(
    data = lombardy_sf,
    fillColor = "white",
    fillOpacity = 0.55,
    color = "#2c3e50",
    weight = 2,
    opacity = 1
  ) %>%
  addCircleMarkers(
    data = stations_sf,
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 6,
    color = "white",
    weight = 1.5,
    fillColor = "#0072B2",
    fillOpacity = 1,
    popup = ~paste0("Station: ", IDStations)
  ) %>%
  leaflet::addControl(
    html = HTML("<b>Agrimonia wind stations in Lombardy</b>"),
    position = "topright"
  )

# ------------------------------------------------------------
# 6. Helper function: save leaflet map as PDF
# ------------------------------------------------------------

save_leaflet_pdf <- function(map, file_stem,
                             width_px = 2400, height_px = 1600,
                             zoom = 2, delay = 2) {
  html_file <- paste0(file_stem, ".html")
  png_file  <- paste0(file_stem, ".png")
  pdf_file  <- paste0(file_stem, ".pdf")
  
  saveWidget(map, html_file, selfcontained = FALSE)
  
  webshot2::webshot(
    url = html_file,
    file = png_file,
    vwidth = width_px,
    vheight = height_px,
    zoom = zoom,
    delay = delay
  )
  
  img <- magick::image_read(png_file)
  magick::image_write(img, path = pdf_file, format = "pdf")
  
  message("Saved: ", pdf_file)
  invisible(pdf_file)
}

# Example export:
# save_leaflet_pdf(
#   map_lombardy,
#   file.path(output_dir, "figures", "lombardy_stations_leaflet"),
#   width_px = 2600,
#   height_px = 1500,
#   zoom = 2,
#   delay = 2
# )

# ------------------------------------------------------------
# 7. Annual station-level mean wind speed for ws10
# ------------------------------------------------------------

wind[, Year := year(Date)]

wind_year <- wind[
  !is.na(ws10_mean),
  .(mean_ws10 = mean(ws10_mean, na.rm = TRUE)),
  by = .(IDStations, Longitude, Latitude, Year)
]

stations_year_sf <- st_as_sf(
  wind_year,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

plot_mean_ws10_year <- ggplot() +
  geom_sf(data = lombardy, fill = "grey98", color = "grey60") +
  geom_sf(data = stations_year_sf, aes(color = mean_ws10), size = 1.8) +
  scale_color_viridis_c(name = "Mean ws10", option = "C") +
  facet_wrap(~Year, ncol = 3) +
  coord_sf() +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 7, face = "bold"),
    plot.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "Annual mean ws10 per station in Lombardy, 2016-2021",
    x = NULL,
    y = NULL
  )

# Example export:
# ggsave(
#   file.path(output_dir, "figures", "annual_mean_per_station_ws10.pdf"),
#   plot_mean_ws10_year,
#   width = 6.0,
#   height = 4.0,
#   device = cairo_pdf
# )

# ------------------------------------------------------------
# 8. Descriptive statistics
# ------------------------------------------------------------

overall_stats <- rbind(
  wind[, .(
    Height = "ws10",
    n = sum(!is.na(ws10_mean)),
    Median = median(ws10_mean, na.rm = TRUE),
    Mean = mean(ws10_mean, na.rm = TRUE),
    IQR = IQR(ws10_mean, na.rm = TRUE),
    SD = sd(ws10_mean, na.rm = TRUE),
    MIN = min(ws10_mean, na.rm = TRUE),
    MAX = max(ws10_mean, na.rm = TRUE)
  )],
  wind[, .(
    Height = "ws100",
    n = sum(!is.na(ws100_mean)),
    Median = median(ws100_mean, na.rm = TRUE),
    Mean = mean(ws100_mean, na.rm = TRUE),
    IQR = IQR(ws100_mean, na.rm = TRUE),
    SD = sd(ws100_mean, na.rm = TRUE),
    MIN = min(ws100_mean, na.rm = TRUE),
    MAX = max(ws100_mean, na.rm = TRUE)
  )]
)

print(overall_stats)

# ------------------------------------------------------------
# 9. Missing-data summary by station
# ------------------------------------------------------------

missing_df <- wind[, .(
  missing_rate_ws10 = mean(is.na(ws10_mean)),
  missing_rate_ws100 = mean(is.na(ws100_mean))
), by = IDStations]

missing_df <- missing_df[order(missing_rate_ws10)]

# ------------------------------------------------------------
# 10. Example time-series plots for selected stations
# ------------------------------------------------------------

example_stations <- c("1264", "STA.IT0591A")

plot_example_ws10 <- ggplot(
  wind[IDStations %in% example_stations],
  aes(x = Date, y = ws10_mean, color = factor(IDStations))
) +
  geom_line(alpha = 0.85, linewidth = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#512DA8")) +
  labs(
    title = "Daily ws10 for selected stations",
    x = "Date",
    y = "ws10 (m/s)",
    color = "Station ID"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10)
  )

# ------------------------------------------------------------
# 11. Helper function: STL deseasonalisation + AR(1) residuals
# ------------------------------------------------------------

prep_stl_ar1 <- function(d, value_col, freq = 365, s_window = "periodic",
                         maxgap = 14, robust = TRUE) {
  
  d <- as.data.table(copy(d))
  d <- d[order(Date)]
  d <- d[, .(Date, y = get(value_col))]
  
  # Build complete daily date grid
  all_dates <- seq(min(d$Date, na.rm = TRUE), max(d$Date, na.rm = TRUE), by = "day")
  full_grid <- data.table(Date = all_dates)
  full_grid <- merge(full_grid, d, by = "Date", all.x = TRUE)
  
  # Interpolate only short gaps
  y_filled <- zoo::na.approx(full_grid$y, x = full_grid$Date, na.rm = FALSE, maxgap = maxgap)
  
  # If gaps remain, return NA columns
  if (anyNA(y_filled)) {
    full_grid[, `:=`(
      seasonal = NA_real_,
      deseason = NA_real_,
      ar1_resid = NA_real_
    )]
    return(full_grid)
  }
  
  # STL decomposition
  y_ts <- ts(y_filled, frequency = freq)
  stl_fit <- tryCatch(
    stl(y_ts, s.window = s_window, robust = robust),
    error = function(e) NULL
  )
  
  if (is.null(stl_fit)) {
    full_grid[, `:=`(
      seasonal = NA_real_,
      deseason = NA_real_,
      ar1_resid = NA_real_
    )]
    return(full_grid)
  }
  
  seasonal_component <- as.numeric(stl_fit$time.series[, "seasonal"])
  deseason_series <- y_filled - seasonal_component
  
  # AR(1) fit on deseasonalised series
  ar1_fit <- tryCatch(
    arima(deseason_series, order = c(1, 0, 0)),
    error = function(e) NULL
  )
  
  if (is.null(ar1_fit)) {
    full_grid[, `:=`(
      seasonal = seasonal_component,
      deseason = deseason_series,
      ar1_resid = NA_real_
    )]
    return(full_grid)
  }
  
  full_grid[, seasonal := seasonal_component]
  full_grid[, deseason := deseason_series]
  full_grid[, ar1_resid := as.numeric(residuals(ar1_fit))]
  
  return(full_grid)
}

# ------------------------------------------------------------
# 12. Helper function: build STL + AR(1) panel for all stations
# ------------------------------------------------------------

build_panel_stl_ar1 <- function(wind_dt, value_col,
                                freq = 365, s_window = "periodic",
                                maxgap = 14, robust = TRUE) {
  
  wind_dt <- as.data.table(copy(wind_dt))
  wind_dt[, Date := as.Date(Date)]
  wind_dt <- wind_dt[!is.na(IDStations) & !is.na(Date)]
  wind_dt[, IDStations := as.character(IDStations)]
  
  coords <- unique(wind_dt[, .(IDStations, Latitude, Longitude)])
  
  out_long <- wind_dt[, {
    prep_stl_ar1(
      .SD,
      value_col = value_col,
      freq = freq,
      s_window = s_window,
      maxgap = maxgap,
      robust = robust
    )
  }, by = IDStations]
  
  out_long <- merge(out_long, coords, by = "IDStations", all.x = TRUE)
  
  raw_wide <- dcast(out_long, Date ~ IDStations, value.var = "y")
  des_wide <- dcast(out_long, Date ~ IDStations, value.var = "deseason")
  res_wide <- dcast(out_long, Date ~ IDStations, value.var = "ar1_resid")
  
  list(
    long = out_long,
    raw_wide = raw_wide,
    des_wide = des_wide,
    res_wide = res_wide,
    raw_mat = as.matrix(raw_wide[, -1]),
    des_mat = as.matrix(des_wide[, -1]),
    res_mat = as.matrix(res_wide[, -1]),
    dates_vec = raw_wide$Date,
    station_ids = colnames(raw_wide)[-1],
    coords = coords
  )
}

# ------------------------------------------------------------
# 13. Helper functions for ACF diagnostics
# ------------------------------------------------------------

acf_long <- function(x, station, lag_max = 60) {
  x <- x[is.finite(x)]
  a <- stats::acf(x, plot = FALSE, lag.max = lag_max)
  
  data.frame(
    Lag = as.numeric(a$lag),
    acf = as.numeric(a$acf),
    Station = station
  ) |>
    dplyr::filter(Lag > 0)
}

plot_acf_facets <- function(df, title_txt, n_ws10, n_ws100) {
  conf_level <- max(1.96 / sqrt(n_ws10), 1.96 / sqrt(n_ws100))
  
  ggplot(df, aes(x = Lag, y = acf, fill = Station)) +
    geom_col(show.legend = FALSE) +
    geom_hline(
      yintercept = c(conf_level, -conf_level),
      linetype = "dashed",
      color = "red"
    ) +
    facet_wrap(~Station, scales = "free_y") +
    scale_fill_manual(values = c("ws10" = "#512DA8", "ws100" = "#00AFBB")) +
    theme_minimal() +
    labs(title = title_txt, x = "Lag", y = "ACF") +
    theme(
      plot.title = element_text(size = 10),
      strip.text = element_text(size = 9),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    )
}

# ------------------------------------------------------------
# 14. Illustrative STL + AR(1) diagnostics for one station
# ------------------------------------------------------------

station_id <- "1266"

station_data <- wind[IDStations == station_id, .(Date, ws10_mean, ws100_mean)]

ws10_obj <- prep_stl_ar1(station_data, "ws10_mean")
ws100_obj <- prep_stl_ar1(station_data, "ws100_mean")

base_theme_ts <- theme_minimal() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10)
  )

plot_raw_ws10 <- ggplot(ws10_obj, aes(Date, y)) +
  geom_line(linewidth = 0.35, color = "#512DA8") +
  labs(title = "ws10 raw", x = NULL, y = "m/s") +
  base_theme_ts

plot_raw_ws100 <- ggplot(ws100_obj, aes(Date, y)) +
  geom_line(linewidth = 0.35, color = "#00AFBB") +
  labs(title = "ws100 raw", x = NULL, y = "m/s") +
  base_theme_ts

plot_deseason_ws10 <- ggplot(ws10_obj, aes(Date, deseason)) +
  geom_line(linewidth = 0.35, color = "#512DA8") +
  labs(title = "ws10 deseasonalised (STL)", x = NULL, y = "m/s") +
  base_theme_ts

plot_deseason_ws100 <- ggplot(ws100_obj, aes(Date, deseason)) +
  geom_line(linewidth = 0.35, color = "#00AFBB") +
  labs(title = "ws100 deseasonalised (STL)", x = NULL, y = "m/s") +
  base_theme_ts

plot_resid_ws10 <- ggplot(ws10_obj, aes(Date, ar1_resid)) +
  geom_line(linewidth = 0.35, color = "#512DA8") +
  labs(title = "ws10 AR(1) residuals", x = NULL, y = NULL) +
  base_theme_ts

plot_resid_ws100 <- ggplot(ws100_obj, aes(Date, ar1_resid)) +
  geom_line(linewidth = 0.35, color = "#00AFBB") +
  labs(title = "ws100 AR(1) residuals", x = NULL, y = NULL) +
  base_theme_ts

figure_stl_ar1_pipeline <-
  (plot_raw_ws10 | plot_raw_ws100) /
  (plot_deseason_ws10 | plot_deseason_ws100) /
  (plot_resid_ws10 | plot_resid_ws100)

# ------------------------------------------------------------
# 15. ACF diagnostic plots for the illustrative station
# ------------------------------------------------------------

lag_max <- 60

n_des_10 <- sum(is.finite(ws10_obj$deseason))
n_des_100 <- sum(is.finite(ws100_obj$deseason))
n_res_10 <- sum(is.finite(ws10_obj$ar1_resid))
n_res_100 <- sum(is.finite(ws100_obj$ar1_resid))

acf_des <- bind_rows(
  acf_long(ws10_obj$deseason, "ws10", lag_max),
  acf_long(ws100_obj$deseason, "ws100", lag_max)
)

acf_des_sq <- bind_rows(
  acf_long(ws10_obj$deseason^2, "ws10", lag_max),
  acf_long(ws100_obj$deseason^2, "ws100", lag_max)
)

acf_res <- bind_rows(
  acf_long(ws10_obj$ar1_resid, "ws10", lag_max),
  acf_long(ws100_obj$ar1_resid, "ws100", lag_max)
)

acf_res_sq <- bind_rows(
  acf_long(ws10_obj$ar1_resid^2, "ws10", lag_max),
  acf_long(ws100_obj$ar1_resid^2, "ws100", lag_max)
)

plot_acf_des <- plot_acf_facets(acf_des, "ACF of deseasonalised series", n_des_10, n_des_100)
plot_acf_des_sq <- plot_acf_facets(acf_des_sq, "ACF of deseasonalised series²", n_des_10, n_des_100)
plot_acf_res <- plot_acf_facets(acf_res, "ACF of AR(1) residuals", n_res_10, n_res_100)
plot_acf_res_sq <- plot_acf_facets(acf_res_sq, "ACF of AR(1) residuals²", n_res_10, n_res_100)

figure_acf_diagnostics <-
  (plot_acf_des | plot_acf_des_sq) /
  (plot_acf_res | plot_acf_res_sq)

# Example export:
# ggsave(
#   file.path(output_dir, "figures", "fig_STL_AR1_pipeline_station1266.pdf"),
#   figure_stl_ar1_pipeline,
#   device = cairo_pdf,
#   width = 7.2,
#   height = 6.0,
#   units = "in"
# )
#
# ggsave(
#   file.path(output_dir, "figures", "fig_ACF_diagnostics_station1266.pdf"),
#   figure_acf_diagnostics,
#   device = cairo_pdf,
#   width = 7.2,
#   height = 4.6,
#   units = "in"
# )

# ------------------------------------------------------------
# 16. Build STL + AR(1) panels for ws10 and ws100
# ------------------------------------------------------------

panel_ws10 <- build_panel_stl_ar1(wind, value_col = "ws10_mean", maxgap = 14)
panel_ws100 <- build_panel_stl_ar1(wind, value_col = "ws100_mean", maxgap = 14)

cat("ws10 residual matrix:", nrow(panel_ws10$res_mat), "x", ncol(panel_ws10$res_mat), "\n")
cat("ws100 residual matrix:", nrow(panel_ws100$res_mat), "x", ncol(panel_ws100$res_mat), "\n")

# ------------------------------------------------------------
# 17. Save processed objects for later modelling steps
# ------------------------------------------------------------

save(
  panel_ws10,
  file = file.path(output_dir, "objects", "wind_STL_AR1_ws10.RData")
)

save(
  panel_ws100,
  file = file.path(output_dir, "objects", "wind_STL_AR1_ws100.RData")
)

cat("Saved STL + AR(1) preprocessing objects for ws10 and ws100.\n")