# ============================================================
# Script: 04_build_spatial_weights.R
# Purpose:
#   Construct spatial weight matrices for Lombardy wind stations:
#   (i) distance-band,
#   (ii) k-nearest neighbours,
#   (iii) directional/advection-based weights.
#
# Inputs:
#   - processed wind data with station coordinates
#   - NUTS shapefile for Lombardy boundary
#
# Outputs:
#   - sparse row-standardised weight matrices
#   - listw objects for Moran's I and diagnostics
#   - optional network plots for the appendix
#
# Author: Ariane Meli Chrisko
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(data.table)
library(sf)
library(spdep)
library(Matrix)
library(ggplot2)

# ------------------------------------------------------------
# 2. Define file paths
# ------------------------------------------------------------

wind_data_file <- "data/processed/wind_prepared.RData"
nuts_file <- "data/raw/NUTS/NUTS_RG_01M_2021_4326.shp"

output_dir <- "results/spatial_weights"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------

load(wind_data_file)   # expected to contain: wind

wind <- as.data.table(wind)
wind[, IDStations := as.character(IDStations)]

nuts2 <- st_read(nuts_file, quiet = TRUE)
lombardy <- nuts2[nuts2$NUTS_ID == "ITC4", ]
lombardy_sf <- st_transform(lombardy, 4326)

# ------------------------------------------------------------
# 4. Prepare unique station locations
# ------------------------------------------------------------

stations <- unique(wind[, .(IDStations, Longitude, Latitude)])
stations_sf <- st_as_sf(stations, coords = c("Longitude", "Latitude"), crs = 4326)
stations_utm <- st_transform(stations_sf, 32632)
coords_utm <- st_coordinates(stations_utm)

# ------------------------------------------------------------
# 5. Helper function: plot neighbour graph
# ------------------------------------------------------------

plot_neighbour_graph <- function(nb_obj, coords_mat, stations_sf_proj, region_sf,
                                 title_text) {
  nb_lines <- spdep::nb2lines(nb_obj, coords = coords_mat)
  nb_lines_sf <- st_as_sf(nb_lines)
  st_crs(nb_lines_sf) <- st_crs(stations_sf_proj)
  
  ggplot() +
    geom_sf(data = region_sf, fill = "grey95", color = "grey60", linewidth = 0.3) +
    geom_sf(data = nb_lines_sf, color = "grey80", linewidth = 0.19, alpha = 0.4) +
    geom_sf(data = stations_sf_proj, color = "black", size = 0.8) +
    ggtitle(title_text) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

# ------------------------------------------------------------
# 6. Distance-band weights
# ------------------------------------------------------------

distance_radius_m <- 55000

nb_distance <- dnearneigh(coords_utm, d1 = 0, d2 = distance_radius_m, longlat = FALSE)

distance_counts <- card(nb_distance)
distance_isolates <- sum(distance_counts == 0)

cat("Distance-band weights\n")
cat("  Min neighbours:", min(distance_counts), "\n")
cat("  Median neighbours:", median(distance_counts), "\n")
cat("  Max neighbours:", max(distance_counts), "\n")
cat("  Isolates:", distance_isolates, "\n")
cat("  Connected components:", n.comp.nb(nb_distance)$nc, "\n")

# Fix isolates by linking each isolated station to its nearest neighbour
if (distance_isolates > 0) {
  isolate_idx <- which(distance_counts == 0)
  nn1 <- knearneigh(coords_utm, k = 1)
  nn_id <- nn1$nn[, 1]
  
  for (i in isolate_idx) {
    nb_distance[[i]] <- nn_id[i]
  }
  
  cat("  After isolate fix:\n")
  cat("  Isolates:", sum(card(nb_distance) == 0), "\n")
  cat("  Connected components:", n.comp.nb(nb_distance)$nc, "\n")
}

lw_distance <- nb2listw(nb_distance, style = "W", zero.policy = TRUE)
W_distance <- Matrix(listw2mat(lw_distance), sparse = TRUE)

distance_plot <- plot_neighbour_graph(
  nb_obj = nb_distance,
  coords_mat = coords_utm,
  stations_sf_proj = stations_utm,
  region_sf = lombardy_sf,
  title_text = "Distance-band neighbours (r = 55 km)"
)

weights_distance_meta <- list(
  W = W_distance,
  type = "distance-band",
  radius_m = distance_radius_m,
  style = "row-standardised",
  crs = "EPSG:32632",
  region = "Lombardy",
  date = Sys.Date()
)

# ------------------------------------------------------------
# 7. k-nearest neighbours weights
# ------------------------------------------------------------

k_neighbours <- 5

nb_knn <- knn2nb(knearneigh(coords_utm, k = k_neighbours, longlat = FALSE))
nb_knn <- make.sym.nb(nb_knn)

knn_counts <- card(nb_knn)

cat("\nk-NN weights\n")
cat("  Min neighbours:", min(knn_counts), "\n")
cat("  Max neighbours:", max(knn_counts), "\n")
cat("  Connected components:", n.comp.nb(nb_knn)$nc, "\n")

lw_knn <- nb2listw(nb_knn, style = "W", zero.policy = TRUE)
W_knn <- Matrix(listw2mat(lw_knn), sparse = TRUE)

knn_plot <- plot_neighbour_graph(
  nb_obj = nb_knn,
  coords_mat = coords_utm,
  stations_sf_proj = stations_utm,
  region_sf = lombardy_sf,
  title_text = "k-nearest neighbours (k = 5)"
)

weights_knn_meta <- list(
  W = W_knn,
  type = "k-nearest neighbours",
  k = k_neighbours,
  style = "row-standardised",
  crs = "EPSG:32632",
  region = "Lombardy",
  date = Sys.Date()
)

# ------------------------------------------------------------
# 8. Directional / advection-based weights
# ------------------------------------------------------------

# Expected direction column in wind data
direction_column <- "WE_mode_wind_direction_100m"

# If the direction column is not already in wind, merge it here from raw data
# or another prepared source before continuing.

stopifnot(direction_column %in% names(wind))

stations_direction <- unique(
  wind[, .(IDStations, Longitude, Latitude, direction_raw = get(direction_column))]
)

compass_to_deg <- function(x) {
  if (is.numeric(x)) return(x %% 360)
  
  x <- toupper(trimws(as.character(x)))
  compass_map <- c(
    "N" = 0, "NNE" = 22.5, "NE" = 45, "ENE" = 67.5,
    "E" = 90, "ESE" = 112.5, "SE" = 135, "SSE" = 157.5,
    "S" = 180, "SSW" = 202.5, "SW" = 225, "WSW" = 247.5,
    "W" = 270, "WNW" = 292.5, "NW" = 315, "NNW" = 337.5
  )
  
  unname(compass_map[x])
}

stations_direction[, dir_deg := vapply(direction_raw, compass_to_deg, numeric(1)) %% 360]
stopifnot(sum(is.na(stations_direction$dir_deg)) == 0)

stations_direction_sf <- st_as_sf(
  stations_direction,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

stations_direction_utm <- st_transform(stations_direction_sf, 32632)
coords_direction_utm <- st_coordinates(stations_direction_utm)

bearing_deg <- function(x1, y1, x2, y2) {
  dx <- x2 - x1
  dy <- y2 - y1
  ang <- atan2(dx, dy) * 180 / pi
  (ang + 360) %% 360
}

cone_half_angle <- 45
max_distance_m <- 100000
distance_decay_scale <- 50000

n_stations <- nrow(coords_direction_utm)
W_directional_dense <- matrix(0, n_stations, n_stations)

for (i in seq_len(n_stations)) {
  xi <- coords_direction_utm[i, 1]
  yi <- coords_direction_utm[i, 2]
  
  for (j in seq_len(n_stations)) {
    if (i == j) next
    
    xj <- coords_direction_utm[j, 1]
    yj <- coords_direction_utm[j, 2]
    
    dij <- sqrt((xi - xj)^2 + (yi - yj)^2)
    if (dij > max_distance_m) next
    
    bij <- bearing_deg(xi, yi, xj, yj)
    ang_diff <- abs(((bij - stations_direction$dir_deg[i] + 180) %% 360) - 180)
    
    if (ang_diff <= cone_half_angle) {
      W_directional_dense[i, j] <- exp(-dij / distance_decay_scale) *
        cos(ang_diff * pi / 180)
    }
  }
}

row_sums_directional <- rowSums(W_directional_dense)
W_directional_rowstd <- W_directional_dense
W_directional_rowstd[row_sums_directional > 0, ] <-
  W_directional_dense[row_sums_directional > 0, ] / row_sums_directional[row_sums_directional > 0]

lw_directional <- mat2listw(W_directional_rowstd, style = "W", zero.policy = TRUE)
W_directional <- Matrix(W_directional_rowstd, sparse = TRUE)

cat("\nDirectional weights\n")
cat("  Density:", mean(W_directional_rowstd > 0), "\n")
cat("  Zero-neighbour stations:", sum(rowSums(W_directional_rowstd) == 0), "\n")

directional_plot <- plot_neighbour_graph(
  nb_obj = lw_directional$neighbours,
  coords_mat = coords_direction_utm,
  stations_sf_proj = stations_direction_utm,
  region_sf = lombardy_sf,
  title_text = sprintf(
    "Directional neighbours (cone = %d°, max = %d km)",
    cone_half_angle,
    round(max_distance_m / 1000)
  )
)

weights_directional_meta <- list(
  W = W_directional,
  type = "directional-advection",
  cone_half_angle = cone_half_angle,
  max_distance_m = max_distance_m,
  distance_decay_scale = distance_decay_scale,
  style = "row-standardised",
  crs = "EPSG:32632",
  region = "Lombardy",
  date = Sys.Date()
)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

saveRDS(weights_distance_meta, file.path(output_dir, "objects", "W_distance_band_55km.rds"))
saveRDS(weights_knn_meta, file.path(output_dir, "objects", "W_knn_k5.rds"))
saveRDS(weights_directional_meta, file.path(output_dir, "objects", "W_directional_advection.rds"))

saveRDS(lw_distance, file.path(output_dir, "objects", "lw_distance_band_55km.rds"))
saveRDS(lw_knn, file.path(output_dir, "objects", "lw_knn_k5.rds"))
saveRDS(lw_directional, file.path(output_dir, "objects", "lw_directional_advection.rds"))

ggsave(
  file.path(output_dir, "figures", "distance_band_graph.pdf"),
  distance_plot,
  width = 6.0,
  height = 4.0,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "figures", "knn_graph.pdf"),
  knn_plot,
  width = 6.0,
  height = 4.0,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "figures", "directional_graph.pdf"),
  directional_plot,
  width = 6.0,
  height = 4.0,
  device = cairo_pdf
)

cat("\nSaved all spatial weight matrices and plots.\n")