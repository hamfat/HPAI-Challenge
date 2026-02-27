rm(list=ls())

theme_pub <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size * 0.9),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size * 0.9),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "bottom"
    )
}

############################################################
## PACKAGES
############################################################

library(tidyverse)
library(sf)
library(lubridate)
library(ggplot2)
library(viridis)
library(patchwork)
library(gganimate)

############################################################
## LOAD DATA
############################################################

population   <- read_csv("population.csv",          show_col_types = FALSE)
cases        <- read_csv("cases.csv",               show_col_types = FALSE)
movement     <- read_csv("movement.csv",            show_col_types = FALSE)
activity     <- read_csv("activity.csv",            show_col_types = FALSE)
prev_culls   <- read_csv("prev_culls.csv",          show_col_types = FALSE)
spatial_risk <- read_csv("spatial_risk_4weeks.csv", show_col_types = FALSE)

############################################################
##  DATE FORMATTING
############################################################

cases <- cases %>%
  mutate(
    date_confirmed  = as.Date(date_confirmed),
    date_suspicious = as.Date(date_suspicious),
    cull_start      = as.Date(cull_start),
    cull_end        = as.Date(cull_end)
  )

movement <- movement %>%
  mutate(date = as.Date(date))

activity <- activity %>%
  mutate(
    date_start = as.Date(date_start),
    date_end   = as.Date(date_end)
  )

prev_culls <- prev_culls %>%
  mutate(
    cull_start = as.Date(cull_start),
    cull_end   = as.Date(cull_end)
  )

############################################################
## LOAD SPATIAL DATA
############################################################

counties  <- st_read("counties_32626.geojson",  quiet = TRUE)
districts <- st_read("districts_32626.geojson", quiet = TRUE)
clc       <- st_read("clc_32626.geojson",       quiet = TRUE)
hrz       <- st_read("hrz_32626.geojson",       quiet = TRUE)

population_sf <- population %>%
  st_as_sf(coords = c("x", "y"), crs = 32626, remove = FALSE)

cases_full <- cases %>%
  left_join(population, by = "farm_id")

cases_sf <- cases_full %>%
  st_as_sf(coords = c("x", "y"), crs = 32626, remove = FALSE)

############################################################
##  TIME FRAME: FIT VS FORECAST
##    Fit: from first confirmed case to 2026-01-13
##    Forecast: next 4 weeks (28 days)
############################################################

start_date      <- min(cases$date_confirmed, na.rm = TRUE)
fit_end_date    <- as.Date("2026-01-13")
forecast_end    <- fit_end_date + 150
sim_dates       <- seq(start_date, forecast_end, by = "day")
n_days          <- length(sim_dates)

############################################################
##  COUNTY DEFINITIONS AND FARM COUNTS
############################################################

county_farms <- population %>%
  count(county, name = "N_c")

counties_vec <- sort(unique(county_farms$county))
C <- length(counties_vec)

county_index <- tibble(
  county = counties_vec,
  idx    = seq_len(C)
)

############################################################
##  COUNTY-LEVEL DAILY CONFIRMED OUTBREAKS
############################################################

cases_county_daily <- cases %>%
  left_join(population %>% select(farm_id, county),
            by = "farm_id") %>%
  filter(!is.na(county)) %>%
  count(date_confirmed, county, name = "n_confirmed")

cases_daily_total <- cases %>%
  count(date_confirmed, name = "n_confirmed") %>%
  rename(date = date_confirmed)

############################################################
##  COUNTY-TO-COUNTY MOVEMENT MATRIX
############################################################

movement_county <- movement %>%
  left_join(population %>% select(farm_id, county),
            by = c("source_farm" = "farm_id")) %>%
  rename(source_county = county) %>%
  left_join(population %>% select(farm_id, county),
            by = c("dest_farm" = "farm_id")) %>%
  rename(dest_county = county) %>%
  filter(!is.na(source_county), !is.na(dest_county)) %>%
  group_by(source_county, dest_county) %>%
  summarise(volume = sum(volume), .groups = "drop")

W_mov <- matrix(
  0,
  nrow = C,
  ncol = C,
  dimnames = list(counties_vec, counties_vec)
)

for (k in seq_len(nrow(movement_county))) {
  sc <- movement_county$source_county[k]
  dc <- movement_county$dest_county[k]
  if (!is.na(sc) && !is.na(dc)) {
    W_mov[sc, dc] <- W_mov[sc, dc] + movement_county$volume[k]
  }
}

row_sums_mov <- rowSums(W_mov)
W_mov[row_sums_mov > 0, ] <- W_mov[row_sums_mov > 0, ] /
  row_sums_mov[row_sums_mov > 0]

############################################################
## COUNTY-TO-COUNTY DISTANCE MATRIX
############################################################

counties_sf <- counties %>%
  filter(county %in% counties_vec) %>%
  arrange(match(county, counties_vec))

county_centroids <- st_centroid(counties_sf)
coords <- st_coordinates(county_centroids)

D <- as.matrix(dist(coords))

dist_scale <- 10000
W_dist <- exp(-D / dist_scale)
diag(W_dist) <- 0

############################################################
## COMBINED TRANSMISSION MATRIX
############################################################

alpha_mov  <- 0.6
alpha_dist <- 0.4

W_comb <- alpha_mov * W_mov + alpha_dist * W_dist

row_sums_comb <- rowSums(W_comb)
W_comb[row_sums_comb > 0, ] <- W_comb[row_sums_comb > 0, ] /
  row_sums_comb[row_sums_comb > 0]

############################################################
##  BUILD COUNTY-DAY CULLING MATRIX
##     Reactive culling delayed by 3 days
##     Preventive culling: completed only, no delay
############################################################

first_confirmed <- min(cases$date_confirmed, na.rm = TRUE)

# Reactive culling only (completed, delayed by 2 days instead of 3)
reactive_culls <- cases %>%
  filter(cull_status == "completed",
         !is.na(cull_start)) %>%
  mutate(
    cull_start_shifted = cull_start + 2   # <-- shifted earlier by 1 day
  ) %>%
  select(farm_id, cull_start = cull_start_shifted)

# Combine (reactive only)
all_culls <- reactive_culls %>%
  mutate(type = "reactive")

# Map farms → counties and to sim day index
all_culls <- all_culls %>%
  left_join(population %>% select(farm_id, county), by = "farm_id") %>%
  filter(!is.na(county)) %>%
  mutate(
    t_idx = match(as.Date(cull_start), sim_dates),
    c_idx = match(county, counties_vec)
  ) %>%
  filter(!is.na(t_idx), !is.na(c_idx))

# Aggregate to county × day
cull_counts <- all_culls %>%
  count(c_idx, t_idx, name = "n_culled")

# Build Cull_mat[county, day]
Cull_mat <- matrix(
  0,
  nrow = C,
  ncol = n_days,
  dimnames = list(counties_vec, as.character(sim_dates))
)

for (k in seq_len(nrow(cull_counts))) {
  i <- cull_counts$c_idx[k]
  t <- cull_counts$t_idx[k]
  Cull_mat[i, t] <- Cull_mat[i, t] + cull_counts$n_culled[k]
}
############################################################
## STOCHASTIC COUNTY-LEVEL SEIR WITH CULLING
############################################################

run_seir_county <- function(beta_within, beta_between, sigma, gamma,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  S_c <- matrix(0, nrow = C, ncol = n_days)
  E_c <- S_c; I_c <- S_c; R_c <- S_c
  
  # Initial susceptible = all farms
  S_c[,1] <- county_farms$N_c[match(counties_vec, county_farms$county)]
  
  # Seeding from observed cases
  for (k in seq_len(nrow(cases))) {
    farm_id_k <- cases$farm_id[k]
    dc        <- cases$date_confirmed[k]
    cty       <- population$county[match(farm_id_k, population$farm_id)]
    if (is.na(cty) || is.na(dc)) next
    
    c_idx <- match(cty, counties_vec)
    t_idx <- match(dc, sim_dates)
    if (is.na(c_idx) || is.na(t_idx)) next
    
    I_c[c_idx, t_idx] <- I_c[c_idx, t_idx] + 1
    if (S_c[c_idx, t_idx] > 0) S_c[c_idx, t_idx] <- S_c[c_idx, t_idx] - 1
  }
  
  # Forward fill
  for (t in 2:n_days) {
    S_c[,t] <- ifelse(S_c[,t] == 0, S_c[,t-1], S_c[,t])
    I_c[,t] <- ifelse(I_c[,t] == 0, I_c[,t-1], I_c[,t])
  }
  
  # Main SEIR loop
  for (t in 2:n_days) {
    N_c <- S_c[,t-1] + E_c[,t-1] + I_c[,t-1] + R_c[,t-1]
    N_c[N_c == 0] <- 1
    
    lambda_within  <- beta_within  * (I_c[,t-1] / N_c)
    lambda_between <- beta_between * as.vector(W_comb %*% (I_c[,t-1] / N_c))
    lambda <- lambda_within + lambda_between
    
    new_E <- rbinom(C, S_c[,t-1], 1 - exp(-lambda))
    new_I <- rbinom(C, E_c[,t-1], sigma)
    new_R <- rbinom(C, I_c[,t-1], gamma)
    
    S_c[,t] <- S_c[,t-1] - new_E
    E_c[,t] <- E_c[,t-1] + new_E - new_I
    I_c[,t] <- I_c[,t-1] + new_I - new_R
    R_c[,t] <- R_c[,t-1] + new_R
    
    # Apply reactive culling (infected farms first)
    to_cull <- Cull_mat[,t]
    
    # Remove from I first
    remove_I <- pmin(I_c[,t], to_cull)
    I_c[,t] <- I_c[,t] - remove_I
    to_cull <- to_cull - remove_I
    
    # If any culling remains, remove from E
   # remove_E <- pmin(E_c[,t], to_cull)
   # E_c[,t] <- E_c[,t] - remove_E
   # to_cull <- to_cull - remove_E
    
    # Susceptible farms should not be culled in reactive-only scenario
    remove_S <- 0
    remove_E <- 0
    
    # Move culled farms to R
    R_c[,t] <- R_c[,t] + remove_I 
    #+ remove_E
  }
  
  list(S_c = S_c, E_c = E_c, I_c = I_c, R_c = R_c)
}

############################################################
##  CALIBRATION (GRID SEARCH) USING STOCHASTIC MODEL
##     Fit to total daily incidence up to fit_end_date
############################################################

obs_daily <- tibble(date = sim_dates) %>%
  left_join(cases_daily_total, by = "date") %>%
  mutate(n_confirmed = replace_na(n_confirmed, 0)) %>%
  filter(date <= fit_end_date)

fit_seir_county <- function(beta_within, beta_between, sigma, gamma) {
  sim <- run_seir_county(beta_within, beta_between, sigma, gamma, seed = 123)
  
  I_c     <- sim$I_c
  I_total <- colSums(I_c)
  I_inc   <- c(0, diff(I_total))
  I_inc[I_inc < 0] <- 0
  
  df <- tibble(
    date     = sim_dates,
    sim_I_inc = I_inc
  ) %>%
    filter(date <= fit_end_date) %>%
    left_join(obs_daily, by = "date")
  
  sum((df$n_confirmed - df$sim_I_inc)^2)
}

beta_within_grid  <- c(0.06, 0.1, 0.2)
beta_between_grid <- c(0.007, 0.1, 0.2)
sigma_grid        <- c(1/3, 1/2)
gamma_grid        <- c(1/7, 1/5)

calib_results <- expand.grid(
  beta_within  = beta_within_grid,
  beta_between = beta_between_grid,
  sigma        = sigma_grid,
  gamma        = gamma_grid
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    ssd = fit_seir_county(beta_within, beta_between, sigma, gamma)
  ) %>%
  ungroup() %>%
  arrange(ssd)

print(calib_results)

best_par <- calib_results %>% slice(1)
print(best_par)

beta_within_best  <- best_par$beta_within
beta_between_best <- best_par$beta_between
sigma_best        <- best_par$sigma
gamma_best        <- best_par$gamma

############################################################
## ENSEMBLE SEIR (MULTI-RUN)
############################################################

run_seir_multi <- function(n_runs,
                           beta_within, beta_between, sigma, gamma) {
  S_arr <- array(0, dim = c(n_runs, C, n_days))
  E_arr <- array(0, dim = c(n_runs, C, n_days))
  I_arr <- array(0, dim = c(n_runs, C, n_days))
  R_arr <- array(0, dim = c(n_runs, C, n_days))
  
  for (r in 1:n_runs) {
    sim <- run_seir_county(beta_within, beta_between, sigma, gamma,
                           seed = 1000 + r)
    S_arr[r,,] <- sim$S_c
    E_arr[r,,] <- sim$E_c
    I_arr[r,,] <- sim$I_c
    R_arr[r,,] <- sim$R_c
  }
  
  list(S = S_arr, E = E_arr, I = I_arr, R = R_arr)
}

n_runs <- 200

multi <- run_seir_multi(
  n_runs       = n_runs,
  beta_within  = beta_within_best,
  beta_between = beta_between_best,
  sigma        = sigma_best,
  gamma        = gamma_best
)

############################################################
## CI (TOTALS ACROSS COUNTIES)
############################################################

get_ci <- function(arr, probs = c(0.025, 0.975)) {
  total <- apply(arr, c(1, 3), sum)  # [run, day]
  tibble(
    date = sim_dates,
    mean = apply(total, 2, mean),
    lo   = apply(total, 2, quantile, probs[1]),
    hi   = apply(total, 2, quantile, probs[2])
  )
}

S_ci <- get_ci(multi$S)
E_ci <- get_ci(multi$E)
I_ci <- get_ci(multi$I)
R_ci <- get_ci(multi$R)

############################################################
## PLOT ENSEMBLE SEIR CI CURVES (WITH FIT/FORECAST SPLIT)
############################################################

plot_ci <- function(df, title, color) {
  ggplot(df, aes(date, mean)) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                fill = color, alpha = 0.25) +
    geom_line(color = color, linewidth = 1.3) +
    geom_vline(xintercept = fit_end_date,
               linetype = "dashed", color = "grey30") +
    annotate(
      "text",
      x = fit_end_date,
      y = max(df$hi, na.rm = TRUE),
      label = "Forecast start",
      hjust = -0.1,
      vjust = 1,
      size = 3
    ) +
    labs(
      title = title,
      x = "Date",
      y = "Number of farms"
    ) +
    theme_pub() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14)
    )
}

pS <- plot_ci(S_ci, "Susceptible (S)", "steelblue")
pE <- plot_ci(E_ci, "Exposed (E)",    "goldenrod")
pI <- plot_ci(I_ci, "Infectious (I)", "firebrick")
pR <- plot_ci(R_ci, "Removed (R)",    "darkgreen")

(pS | pE) / (pI | pR)

ggsave(
  "SEIR_Q6_noPrev.pdf",
  plot = (pS | pE) / (pI | pR),
  width = 12,
  height = 10,
  device = cairo_pdf
)

############################################################
##  SPECIES & PRODUCTION COMPOSITION WITHIN COUNTIES
############################################################

species_county <- population %>%
  filter(!is.na(county), !is.na(species)) %>%
  count(county, species, name = "n_farms") %>%
  group_by(county) %>%
  mutate(
    total_farms = sum(n_farms),
    prop        = n_farms / total_farms
  ) %>%
  ungroup()

production_county <- population %>%
  filter(!is.na(county), !is.na(production)) %>%
  count(county, production, name = "n_farms") %>%
  group_by(county) %>%
  mutate(
    total_farms = sum(n_farms),
    prop        = n_farms / total_farms
  ) %>%
  ungroup()

############################################################
## SPECIES-SPECIFIC ENSEMBLE SEIR (I ONLY)
############################################################

species_levels <- unique(species_county$species)

species_seir_ci <- function(arr, species_name) {
  w <- species_county %>%
    filter(species == species_name) %>%
    right_join(tibble(county = counties_vec), by = "county") %>%
    arrange(county) %>%
    pull(prop)
  w[is.na(w)] <- 0
  
  total <- sapply(1:dim(arr)[1], function(r) {
    apply(arr[r,,], 2, function(x) sum(x * w))
  })
  total <- t(total)  # [run, day]
  
  tibble(
    date    = sim_dates,
    mean    = apply(total, 2, mean),
    lo      = apply(total, 2, quantile, 0.025),
    hi      = apply(total, 2, quantile, 0.975),
    species = species_name
  )
}

I_species_ci <- map_dfr(species_levels, ~ species_seir_ci(multi$I, .x))

species_cols <- c(
  "chicken" = "#1f78b4",
  "duck"    = "#e31a1c"
)

p_I_species <- ggplot(
  I_species_ci,
  aes(date, mean, color = species, fill = species)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_vline(xintercept = fit_end_date,
             linetype = "dashed", color = "grey30") +
  scale_color_manual(values = species_cols) +
  scale_fill_manual(values = species_cols) +
  labs(
    title = "Species-Specific Infectious Farms",
    x = "Date",
    y = "Number of infectious farms",
    color = "Species",
    fill  = "Species"
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_I_species)

ggsave(
  "Ispecies_Q5.pdf",
  plot = p_I_species,
  width = 12,
  height = 8,
  device = cairo_pdf
)

############################################################
## PRODUCTION-SPECIFIC ENSEMBLE SEIR (I ONLY)
############################################################

prod_levels <- unique(production_county$production)

prod_seir_ci <- function(arr, prod_name) {
  w <- production_county %>%
    filter(production == prod_name) %>%
    right_join(tibble(county = counties_vec), by = "county") %>%
    arrange(county) %>%
    pull(prop)
  w[is.na(w)] <- 0
  
  total <- sapply(1:dim(arr)[1], function(r) {
    apply(arr[r,,], 2, function(x) sum(x * w))
  })
  total <- t(total)
  
  tibble(
    date      = sim_dates,
    mean      = apply(total, 2, mean),
    lo        = apply(total, 2, quantile, 0.025),
    hi        = apply(total, 2, quantile, 0.975),
    production = prod_name
  )
}

I_prod_ci <- map_dfr(prod_levels, ~ prod_seir_ci(multi$I, .x))

production_cols <- c(
  "broiler_1"    = "#1f78b4",
  "broiler_2"    = "#33a02c",
  "conventional" = "#e31a1c",
  "layer"        = "#ff7f00",
  "organic"      = "#6a3d9a"
)

p_I_prod <- ggplot(
  I_prod_ci,
  aes(date, mean, color = production, fill = production)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_vline(xintercept = fit_end_date,
             linetype = "dashed", color = "grey30") +
  scale_color_manual(values = production_cols) +
  scale_fill_manual(values = production_cols) +
  labs(
    title = "Production-Type Specific Infectious Farms",
    x = "Date",
    y = "Number of infectious farms",
    color = "Production",
    fill  = "Production"
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_I_prod)

ggsave(
  "Iprod_Q6_noPrev.pdf",
  plot = p_I_prod,
  width = 12,
  height = 8,
  device = cairo_pdf
)

############################################################
##  COUNTY-LEVEL ENSEMBLE I(t): MEAN + CI
############################################################

get_county_I_ci <- function(arr) {
  out <- vector("list", C)
  for (j in 1:C) {
    total <- arr[, j, ]  # [run, day]
    out[[j]] <- tibble(
      date   = sim_dates,
      mean   = apply(total, 2, mean),
      lo     = apply(total, 2, quantile, 0.025),
      hi     = apply(total, 2, quantile, 0.975),
      county = counties_vec[j]
    )
  }
  bind_rows(out)
}

I_county_ci <- get_county_I_ci(multi$I)

Icounty<-ggplot(I_county_ci, aes(date, mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.15) +
  geom_line(color = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = fit_end_date,
             linetype = "dashed", color = "grey30") +
  facet_wrap(~ county, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "County-level infectious farms",
    x = "Date",
    y = "Number of infectious farms"
  )+theme_pub()

print(Icounty)

ggsave(
  "Icounty_Q6_noPrev.pdf",
  plot = Icounty,
  width = 12,
  height = 8,
  device = cairo_pdf
)
############################################################
## JOIN ENSEMBLE I(t) WITH COUNTY POLYGONS (MAPS + ANIMATION)
############################################################

I_county_map <- counties_sf %>%
  left_join(I_county_ci, by = "county")

counties_centroids <- counties_sf %>%
  st_centroid()

day_to_plot <- 5

p_single_day <- ggplot(
  I_county_map %>% filter(date == sim_dates[day_to_plot])
) +
  geom_sf(aes(fill = mean), color = "white", size = 0.2) +
  geom_sf_text(
    data = counties_centroids,
    aes(label = county),
    size = 3,
    color = "white"
  ) +
  scale_fill_viridis_c(option = "C") +
  theme_pub() +
  theme(
    legend.position = "right"
  ) +
  labs(
    title = paste("Mean infectious farms -", sim_dates[day_to_plot]),
    fill  = "I(t)",
    x = NULL,
    y = NULL
  )

print(p_single_day)

ggsave(
  "infectiouscountyD1.pdf",
  plot = p_single_day,
  width = 12,
  height = 10,
  device = cairo_pdf
)

days_to_plot <- sim_dates[seq(1, n_days, by = 10)]

p_facet2 <- ggplot(
  I_county_map %>% filter(date %in% days_to_plot)
) +
  geom_sf(aes(fill = mean), color = "white") +
  geom_sf_text(
    data = counties_centroids,
    aes(label = county),
    size = 2.8,
    color = "white",
    check_overlap = TRUE
  ) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ date) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(
    title = "Mean infectious farms over time",
    fill  = "I(t)"
  )

print(p_facet2)

ggsave(
  "infectiouscountyovertime_Q5.pdf",
  plot = p_facet2,
  width = 12,
  height = 10,
  device = cairo_pdf
)


early_dates <- days_to_plot[1:3]

p_early <- ggplot(
  I_county_map %>% filter(date %in% early_dates)
) +
  geom_sf(aes(fill = mean), color = "white") +
  geom_sf_text(
    data = counties_centroids,
    aes(label = county),
    size = 2.8,
    color = "white",
    check_overlap = TRUE
  ) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ date) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(
    title = "Mean infectious farms — early outbreak period",
    fill  = "I(t)"
  )
print(p_early)



early_dates  <- days_to_plot[1:3]
later_dates  <- days_to_plot[4:6]

# Early-date plot with its own color scale
p_early <- ggplot(
  I_county_map %>% filter(date %in% early_dates)
) +
  geom_sf(aes(fill = mean), color = "white") +
  geom_sf_text(
    data = counties_centroids,
    aes(label = county),
    size = 2.8,
    color = "white",
    check_overlap = TRUE
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt") +
  facet_wrap(~ date) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(
    title = "Early outbreak period",
    fill  = "I(t)"
  )

# Later-date plot with original scale
p_late <- ggplot(
  I_county_map %>% filter(date %in% later_dates)
) +
  geom_sf(aes(fill = mean), color = "white") +
  geom_sf_text(
    data = counties_centroids,
    aes(label = county),
    size = 2.8,
    color = "white",
    check_overlap = TRUE
  ) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ date) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(
    title = "Later outbreak period",
    fill  = "I(t)"
  )

# Combine into one 6‑panel figure
p_all <- p_early / p_late

print(p_all)

ggsave(
  filename = "infectious_maps_Q5.pdf",
  plot = p_all,
  device = cairo_pdf,
  width = 12,
  height = 10,
  units = "in"
)


p_anim <- ggplot(I_county_map) +
  geom_sf(aes(fill = mean), color = NA) +
  scale_fill_viridis_c(option = "C") +
  theme_minimal() +
  labs(
    title = "Ensemble infectious farms - {frame_time}",
    fill  = "I(t)"
  ) +
  transition_time(date) +
  ease_aes("linear")

anim <- animate(
  p_anim,
  fps = 10,
  width = 900,
  height = 700,
  renderer = gifski_renderer()
)

anim_save("ensemble_I_animation.gif", animation = anim)

anim <- animate(
  p_anim,
  fps = 10,
  width = 900,
  height = 700,
  renderer = av_renderer()
)

anim_save("ensemble_I_animation.mp4", animation = anim)

############################################################
## ENSEMBLE TOTAL INCIDENCE VS OBSERVED (WITH FIT/FORECAST SPLIT)
############################################################

I_total <- apply(multi$I, c(1, 3), sum)

I_incidence_runs <- t(apply(I_total, 1, function(x) {
  inc <- c(0, diff(x))
  inc[inc < 0] <- 0
  inc
}))

I_incidence_ci <- tibble(
  date = sim_dates,
  mean = apply(I_incidence_runs, 2, mean),
  lo   = apply(I_incidence_runs, 2, quantile, 0.025),
  hi   = apply(I_incidence_runs, 2, quantile, 0.975)
)

obs_incidence <- cases %>%
  count(date_confirmed, name = "obs_inc") %>%
  rename(date = date_confirmed) %>%
  mutate(obs_inc = replace_na(obs_inc, 0))

compare_inc <- I_incidence_ci %>%
  left_join(obs_incidence, by = "date") %>%
  mutate(obs_inc = replace_na(obs_inc, 0))

p_compare_inc <- ggplot(compare_inc, aes(date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = mean),
            color = "steelblue", linewidth = 1.3) +
  geom_point(aes(y = if_else(date <= fit_end_date, obs_inc, NA_real_)),
             color = "black", size = 2) +
  geom_line(aes(y = if_else(date <= fit_end_date, obs_inc, NA_real_)),
            color = "black", linewidth = 1) +
  geom_vline(xintercept = fit_end_date,
             linetype = "dashed", color = "grey30") +
  labs(
    title = "Observed Incidence vs Predicted Incidence",
    #subtitle = "Daily new infectious farms: ensemble mean + 95% CI vs confirmed cases",
    x = "Date",
    y = expression("New infectious farms per day")
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_compare_inc)

ggsave(
  "predictedvsobserved_incidence_Q5.pdf",
  plot = p_compare_inc,
  width = 10,
  height = 6,
  device = cairo_pdf
)

############################################################
##  ACTUAL OBSERVED INFECTIOUS FARMS PER COUNTY PER DAY
############################################################

actual_I <- cases %>%
  left_join(population %>% select(farm_id, county), by = "farm_id") %>%
  count(date_confirmed, county, name = "I_actual") %>%
  rename(date = date_confirmed)

actual_I_full <- expand_grid(
  date = sim_dates,
  county = counties_vec
) %>%
  left_join(actual_I, by = c("date", "county")) %>%
  mutate(I_actual = replace_na(I_actual, 0))

actual_I_map <- counties_sf %>%
  left_join(actual_I_full, by = "county")

day_to_plot <- 6

p_actual_single <- ggplot(
  actual_I_map %>% filter(date == sim_dates[day_to_plot])
) +
  geom_sf(aes(fill = I_actual), color = "grey40", size = 0.2) +
  scale_fill_viridis_c(option = "C") +
  theme_minimal() +
  labs(
    title = paste("Actual infectious farms -", sim_dates[day_to_plot]),
    fill  = "I_actual"
  )

print(p_actual_single)

days_to_plot <- sim_dates[seq(1, n_days, by = 10)]

p_actual_facet <- ggplot(
  actual_I_map %>% filter(date %in% days_to_plot)
) +
  geom_sf(aes(fill = I_actual), color = NA) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ date) +
  theme_minimal() +
  labs(
    title = "Actual infectious farms over time",
    fill  = "I_actual"
  )

print(p_actual_facet)

p_actual_anim <- ggplot(actual_I_map) +
  geom_sf(aes(fill = I_actual), color = NA) +
  scale_fill_viridis_c(option = "C") +
  theme_minimal() +
  labs(
    title = "Actual infectious farms - {frame_time}",
    fill  = "I_actual"
  ) +
  transition_time(date) +
  ease_aes("linear")

anim <- animate(
  p_actual_anim,
  fps = 5,
  width = 900,
  height = 700,
  renderer = gifski_renderer()
)

anim_save("actual_I_animation.gif", animation = anim)

############################################################
##  COUNTY-LEVEL ACTUAL VS ENSEMBLE I(t) WITH CI
############################################################

get_county_I_mean <- function(arr) {
  out <- vector("list", C)
  for (j in 1:C) {
    total <- arr[, j, ]  # [run, day]
    out[[j]] <- tibble(
      date      = sim_dates,
      I_ensemble = apply(total, 2, mean),
      county    = counties_vec[j]
    )
  }
  bind_rows(out)
}

I_county_mean <- get_county_I_mean(multi$I)

county_ts_compare <- actual_I_full %>%
  left_join(I_county_mean, by = c("county", "date")) %>%
  select(county, date, I_actual, I_ensemble)

I_county_ci_band <- get_county_I_ci(multi$I) %>%
  select(date, lo, hi, county)

county_ts_compare_ci <- county_ts_compare %>%
  left_join(I_county_ci_band, by = c("county", "date"))

p_county_compare <- ggplot(county_ts_compare_ci, aes(date)) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    fill = "firebrick",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = I_actual),
    color = "black",
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = I_ensemble),
    color = "firebrick",
    linewidth = 0.9
  ) +
  geom_vline(xintercept = fit_end_date,
             linetype = "dashed", color = "grey30") +
  facet_wrap(~ county, scales = "free_y") +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "County-level Infectious Farms: Actual vs Predicted (95% CI)",
    x = "Date",
    y = "Number of infectious farms",
    caption = "Black = Actual | Red = Predicted mean | Ribbon = 95% CI"
  )

print(p_county_compare)

ggsave(
  "Icountyactualvensemble_Q5.pdf",
  plot = p_county_compare,
  width = 12,
  height = 10,
  device = cairo_pdf
)

############################################################
## PREDICTED BURDEN OVER 4-WEEK HORIZON (FROM FIT END)
############################################################

pred_start <- fit_end_date + 1
pred_end   <- fit_end_date + 28

pred_inc_county <- tibble(
  date   = rep(sim_dates, each = C),
  county = rep(counties_vec, times = n_days),
  I      = as.vector(apply(multi$I, c(2,3), mean))  # ensemble mean I
) %>%
  group_by(county) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(I_prev = lag(I, default = first(I))) %>%
  ungroup() %>%
  mutate(I_inc = pmax(I - I_prev, 0)) %>%
  filter(date >= pred_start, date <= pred_end) %>%
  group_by(county) %>%
  summarise(
    total_new_I = sum(I_inc),
    .groups = "drop"
  )

print(pred_inc_county)

counties_pred <- counties_sf %>%
  left_join(pred_inc_county, by = "county") %>%
  mutate(total_new_I = replace_na(total_new_I, 0))

counties_pred_centroids <- counties_pred %>%
  st_centroid()

p_pred_map <- ggplot(counties_pred) +
  geom_sf(aes(fill = total_new_I), color = "white", size = 0.2) +
  geom_sf_text(
    data = counties_pred_centroids,
    aes(label = county),
    size = 3.2,
    color = "white"
  ) +
  scale_fill_viridis_c(
    option = "C",
    name = expression("New infections")
  ) +
  labs(
    title = "Predicted New Infectious Farms Over Next 4 Weeks",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

print(p_pred_map)

ggsave(
  "predinfecmap_Q5.png",
  plot = p_pred_map,
  width = 10,
  height = 7,
  dpi = 600
)

ggsave(
  "predinfecmap_Q5.pdf",
  plot = p_pred_map,
  width = 10,
  height = 7,
  device = cairo_pdf
)

############################################################
##  SPECIES- AND PRODUCTION-WEIGHTED RISK MAPS
############################################################

species_risk <- species_county %>%
  left_join(pred_inc_county, by = "county") %>%
  mutate(
    total_new_I       = replace_na(total_new_I, 0),
    exp_new_I_species = total_new_I * prop
  )

top_risk_counties <- species_risk %>%
  group_by(county) %>%
  summarise(total_new_I = max(total_new_I), .groups = "drop") %>%
  arrange(desc(total_new_I)) %>%
  slice_head(n = 14) %>%
  pull(county)

species_risk_top <- species_risk %>%
  filter(county %in% top_risk_counties) %>%
  arrange(desc(total_new_I), county, species)

print(
  species_risk_top %>%
    select(county, species, total_new_I, exp_new_I_species)
)

species_cols <- c(
  "chicken" = "#1f78b4",
  "duck"    = "#e31a1c"
)

p_species_risk <- ggplot(
  species_risk %>% filter(county %in% top_risk_counties),
  aes(x = species, y = exp_new_I_species, fill = species)
) +
  geom_col(width = 0.75, color = "white") +
  facet_wrap(~ county, scales = "free_y") +
  scale_fill_manual(values = species_cols) +
  labs(
    title = "Species-Specific Risk by County",
    x = "Species",
    y = expression("Expected new infectious farms"),
    fill = "Species"
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_species_risk)

ggsave(
  "speciesrisk_Q5.pdf",
  plot = p_species_risk,
  width = 10,
  height = 7,
  device = cairo_pdf
)

production_risk <- production_county %>%
  left_join(pred_inc_county, by = "county") %>%
  mutate(
    total_new_I    = replace_na(total_new_I, 0),
    exp_new_I_prod = total_new_I * prop
  )

production_risk_top <- production_risk %>%
  filter(county %in% top_risk_counties) %>%
  arrange(desc(total_new_I), county, production)

print(
  production_risk_top %>%
    select(county, production, total_new_I, exp_new_I_prod)
)

production_cols <- c(
  "broiler_1"    = "#1f78b4",
  "broiler_2"    = "#33a02c",
  "conventional" = "#e31a1c",
  "layer"        = "#ff7f00",
  "organic"      = "#6a3d9a"
)

p_production_risk <- ggplot(
  production_risk %>% filter(county %in% top_risk_counties),
  aes(x = production, y = exp_new_I_prod, fill = production)
) +
  geom_col(width = 0.75, color = "white") +
  facet_wrap(~ county, scales = "free_y") +
  scale_fill_manual(values = production_cols) +
  labs(
    title = "Production-Type Specific Risk by County",
    x = "Production type",
    y = expression("Expected new infectious farms (" * Delta * I * ")"),
    fill = "Production"
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_production_risk)

ggsave(
  "prodrisk_Q5.pdf",
  plot = p_production_risk,
  width = 10,
  height = 7,
  device = cairo_pdf
)

dominant_species <- species_county %>%
  group_by(county) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(county, dom_species = species, dom_species_prop = prop)

dominant_production <- production_county %>%
  group_by(county) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(county, dom_prod = production, dom_prod_prop = prop)

counties_rich <- counties_sf %>%
  left_join(pred_inc_county,     by = "county") %>%
  left_join(dominant_species,    by = "county") %>%
  left_join(dominant_production, by = "county") %>%
  mutate(
    total_new_I      = replace_na(total_new_I, 0),
    dom_species_prop = replace_na(dom_species_prop, 0),
    dom_prod_prop    = replace_na(dom_prod_prop, 0)
  )

p_dom_species <- ggplot(counties_rich) +
  geom_sf(aes(fill = total_new_I), color = "white", size = 0.2) +
  geom_sf_text(aes(label = dom_species),
               size = 3.2, color = "white") +
  scale_fill_viridis_c(
    option = "C",
    name = expression("New infections")
  ) +
  labs(
    title = "Predicted County Burden by Dominant Species",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

print(p_dom_species)

ggsave(
  "dom-species_Q5.pdf",
  plot = p_dom_species,
  width = 10,
  height = 8,
  device = cairo_pdf
)

p_dom_prod <- ggplot(counties_rich) +
  geom_sf(aes(fill = total_new_I), color = "white", size = 0.2) +
  geom_sf_text(aes(label = dom_prod),
               size = 3.2, color = "white") +
  scale_fill_viridis_c(
    option = "C",
    name = expression("New infections")
  ) +
  labs(
    title = "Predicted County Burden by Dominant Production-Type",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

print(p_dom_prod)

ggsave(
  "dom-prod_5.pdf",
  plot = p_dom_prod,
  width = 10,
  height = 8,
  device = cairo_pdf
)


