library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(sf)
library(FNN)
library(purrr)


#### plotting function

theme_epi <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Titles
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0),
      
      # Axes
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text  = element_text(size = base_size),
      
      # Gridlines
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10),
      
      # Legend
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text  = element_text(size = base_size)
    )
}

############################################################
# Load data 
############################################################

population <- read_csv("population.csv", show_col_types = FALSE)
cases      <- read_csv("cases.csv",      show_col_types = FALSE)
movement   <- read_csv("movement.csv",   show_col_types = FALSE)
activity   <- read_csv("activity.csv",   show_col_types = FALSE)
prev_culls <- read_csv("prev_culls.csv", show_col_types = FALSE)
spatial_risk <- read_csv("spatial_risk_4weeks.csv", show_col_types = FALSE)

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
# Load spatial data
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
##  Descriptive epidemiology
############################################################

tab_species_prod <- cases_full %>%
  count(species, production, name = "n_outbreaks")
print(tab_species_prod)

daily_incidence_obs <- cases %>%
  count(date_confirmed, name = "n_confirmed")

p_daily_obs <- ggplot(daily_incidence_obs, aes(date_confirmed, n_confirmed)) +
  geom_col(fill = "red") +
  labs(
    title = "Observed Daily Confirmed Outbreaks",
    x = "Date",
    y = "New Outbreaks"
  ) +
  theme_epi()

print(p_daily_obs)

ggsave(
  "EpiCurve_DailyCases.png",
  plot = p_daily_obs,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)

weekly_incidence_obs <- cases %>%
  mutate(week = floor_date(date_confirmed, "week")) %>%
  count(week, name = "n_confirmed")

p_weekly_obs <- ggplot(weekly_incidence_obs, aes(week, n_confirmed)) +
  geom_col(fill = "blue") +
  theme_minimal() +
  labs(title = "Observed weekly confirmed outbreaks",
       x = "Week", y = "New outbreaks")+
  theme_epi()
print(p_weekly_obs)


ggsave(
  "EpiCurve_WeeklyCases.png",
  plot = p_weekly_obs,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)

cases_species <- cases_full %>%
  count(species, name = "n_outbreaks")

cases_specie<-ggplot(cases_species, aes(x = species, y = n_outbreaks, fill = species)) +
  geom_col() +
  scale_fill_manual(values = c(chicken = "darkorange", duck = "dodgerblue")) +
  theme_minimal() +
  labs(
    title = "Total Outbreaks by Species",
    x = "Species",
    y = "Number of Outbreaks"
  )+theme_epi()

print(cases_specie)

ggsave(
  "EpiCurve_Species.png",
  plot = cases_specie,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)

cases_prod <- cases_full %>%
  count(production, name = "n_outbreaks")

case_prodT <- ggplot(cases_prod, aes(
  x = reorder(production, n_outbreaks),
  y = n_outbreaks,
  fill = production
)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Total Outbreaks by Production Type",
    x = "Production Type",
    y = "Number of Outbreaks"
  ) +
  theme_epi()
print(case_prodT)

ggsave(
  "EpiCurve_ProdType.png",
  plot = case_prodT,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)


cases_species_time <- cases_full %>%
  count(date_confirmed, species, name = "n_outbreaks")

cases_species_time1<-ggplot(cases_species_time, aes(x = date_confirmed, y = n_outbreaks, fill = species)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(chicken = "darkorange", duck = "dodgerblue")) +
  theme_minimal() +
  labs(
    title = "Temporal Evolution of Outbreaks by Species",
    x = "Date",
    y = "New Outbreaks"
  )+
  theme_epi()
print(cases_species_time1)

ggsave(
  "EpiCurve_Species1.png",
  plot = cases_species_time1,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)


cases_prod_time <- cases_full %>%
  count(date_confirmed, production, name = "n_outbreaks")

cases_prodT1<-ggplot(cases_prod_time, aes(x = date_confirmed, y = n_outbreaks, fill = production)) +
  geom_col() +
  facet_wrap(~ production, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Temporal Evolution of Outbreaks by Production Type",
    x = "Date",
    y = "New Outbreaks"
  )+
  theme_epi()
print(cases_prodT1)

ggsave(
  "EpiCurve_prodT1.png",
  plot = cases_prodT1,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)



cases_by_district <- cases_full %>%
  count(district, name = "n_outbreaks")

districts_cases <- districts %>%
  left_join(cases_by_district, by = "district") %>%
  mutate(n_outbreaks = replace_na(n_outbreaks, 0))

p_map_points <- ggplot() +
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  geom_sf(data = hrz, fill = NA, color = "red") +
  geom_sf(data = cases_sf, aes(color = species), size = 2) +
  scale_color_manual(values = c(chicken = "darkorange", duck = "dodgerblue")) +
  theme_minimal() +
  labs(title = "Spatial distribution of outbreaks",
       color = "Species")+
  theme_epi()

print(p_map_points)


ggsave(
  "casesSpatial.png",
  plot = p_map_points,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)



spatial_ProdT<-ggplot() +
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  geom_sf(data = hrz, fill = NA, color = "red") +
  geom_sf(data = cases_sf, aes(color = production), size = 2) +
  theme_minimal() +
  labs(
    title = "Spatial Distribution of Outbreaks by Production Type",
    color = "Production Type"
  )+
  theme_epi()

print(spatial_ProdT)


ggsave(
  "casesSpatialProdT.png",
  plot = spatial_ProdT,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)



p_map_choro <- ggplot(districts_cases) +
  geom_sf(aes(fill = n_outbreaks), color = NA) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Observed outbreaks per district",
       fill = "Outbreaks")+theme_epi()
print(p_map_choro)

ggsave(
  "casesSpatialDistrict.png",
  plot = p_map_choro,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)


cases_by_district_prod <- cases_full %>%
  count(district, production, name = "n_outbreaks")

districts_prod <- districts %>%
  left_join(cases_by_district_prod, by = "district") %>%
  mutate(n_outbreaks = replace_na(n_outbreaks, 0))

ggplot() +
  # full island background
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  
  # district-level outbreak shading
  geom_sf(data = districts_prod, aes(fill = n_outbreaks), color = NA) +
  
  # facet by production type
  facet_wrap(~ production) +
  
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Outbreaks per District by Production Type (Full Island Map)",
    fill = "Outbreaks"
  )

cases_by_district_prod <- cases_full %>%
  count(district, production, name = "n_outbreaks")

districts_prod <- districts %>%
  left_join(cases_by_district_prod, by = "district") %>%
  mutate(n_outbreaks = replace_na(n_outbreaks, 0))

SpatialProdDist<-ggplot() +
  # full island background
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  
  # HRZ boundaries
  geom_sf(data = hrz, fill = NA, color = "red", linewidth = 0.7) +
  
  # district-level outbreak shading
  geom_sf(data = districts_prod, aes(fill = n_outbreaks), color = NA, alpha = 0.9) +
  
  # facet by production type
  facet_wrap(~ production) +
  
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Outbreaks per District by Production Type",
    fill = "Outbreaks"
  )+theme_epi()
print(SpatialProdDist)

ggsave(
  "SpatialProdDist.png",
  plot = SpatialProdDist,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)



cases_by_district_prod <- cases_full %>%
  count(district, production, name = "n_outbreaks")

districts_prod <- districts %>%
  left_join(cases_by_district_prod, by = "district") %>%
  mutate(n_outbreaks = replace_na(n_outbreaks, 0))

SpatialSpeciesDist<-ggplot() +
  # Full island background
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  
  # HRZ boundaries
  geom_sf(data = hrz, fill = NA, color = "red", linewidth = 0.7) +
  
  # District shading by outbreak count
  geom_sf(data = districts_prod, aes(fill = n_outbreaks), color = NA, alpha = 0.9) +
  
  # Outbreak points
  geom_sf(data = cases_sf, aes(color = species), size = 1.8, alpha = 0.9) +
  
  # Facet by production type
  facet_wrap(~ production) +
  
  scale_fill_viridis_c() +
  scale_color_manual(values = c(chicken = "darkorange", duck = "dodgerblue")) +
  
  theme_minimal() +
  labs(
    title = "Outbreaks per District by Prod Type by Speices",
    fill = "Outbreaks",
    color = "Species"
  )+theme_epi()
print(SpatialSpeciesDist)

ggsave(
  "SpatialSpecieProD.png",
  plot = SpatialSpeciesDist,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)

county_centroids <- st_centroid(counties)


#district_centroids <- st_centroid(districts)

Spatialcountycases<-ggplot() +
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  geom_sf(data = hrz, fill = NA, color = "yellow", linewidth = 0.7) +
  geom_sf(data = districts_prod, aes(fill = n_outbreaks), color = NA, alpha = 0.9) +
  geom_sf(data = cases_sf, aes(color = species), size = 1.8, alpha = 0.9) +
  
  # District labels
  #geom_sf_text(data = district_centroids, aes(label = district), size = 3, color = "black") +
  #county
  geom_sf_text(data = county_centroids, aes(label = county), size = 4, color = "white")

facet_wrap(~ production) +
  
  scale_fill_viridis_c() +
  scale_color_manual(values = c(chicken = "darkorange", duck = "dodgerblue")) +
  
  theme_minimal() +
  labs(
    title = "Outbreaks per District by Production Type (with District Labels)",
    fill = "Outbreaks",
    color = "Species"
  )+theme_epi()
print(Spatialcountycases)



ggsave(
  "Spatialcountycases.png",
  plot = Spatialcountycases ,
  dpi = 600,
  width = 12,
  height = 8,
  units = "in"
)

##########################################
#Map of 1 km risk farms
###########################################
population_sf_risk <- population_sf %>%
  left_join(spatial_risk %>% select(farm_id, within_1km),
            by = "farm_id")

p_risk_map <- ggplot() +
  geom_sf(data = counties, fill = "grey95", color = "grey70") +
  geom_sf(data = hrz, fill = NA, color = "red") +
  geom_sf(data = population_sf_risk,
          aes(color = within_1km), size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "blue")) +
  theme_minimal() +
  labs(title = "Farms within 1 km of infected farms",
       color = "Within 1 km")
print(p_risk_map)




##############################
##Spatial temporal evolution
##############################
library(cowplot)

counties = read_sf("counties_32626.geojson")
districts = read_sf("districts_32626.geojson")
ggplot(data = counties) +
  geom_sf(fill = "blue", color = "white") +
  theme_minimal()


population = read.csv("population.csv",header=T)
head(population)
population_sf = st_as_sf(x = population, coords = c("x", "y"), crs = "EPSG:32626")
cases = read.csv("cases.csv",header=T)
head(cases)
tail(cases)
cases %>% mutate(date_suspicious = ifelse(date_suspicious=="",date_confirmed,date_suspicious),
                 day_sus = as.numeric(as.Date(date_suspicious)-as.Date("2025-12-19")),
                 day_conf = as.numeric(as.Date(date_confirmed)-as.Date("2025-12-19")),
                 cull_end = ifelse(cull_end=="","2026-01-18",cull_end),
                 day_cull_end = as.numeric(as.Date(cull_end)-as.Date("2025-12-19"))) -> cases

cases_daily = read.csv("HPAI_farms_infection_start_end day.csv",header=T)
head(cases_daily)
cases_daily %>% select(c(farm_id:Day25)) %>% 
  mutate(across(Day1:Day26, ~NA)) -> cases_daily

cases_daily %>% pivot_longer(cols = starts_with("Day"),      # Selects columns starting with "Year"
                             names_to = "Day",               # New column for original column names
                             values_to = "Present" )%>%
  mutate(Day_numeric = as.numeric(str_replace(Day,"Day","")))-> cases_daily_long

cases_daily_long$Sus = NA

cases_daily_long %>% left_join(cases,by="farm_id") -> cases_daily_long

cases_daily_long %>%  rowwise() %>% mutate(
  Present = case_when(is.na(day_conf)==F && is.na(day_sus)==F &&
                        Day_numeric >=day_sus &&
                        Day_numeric <day_conf ~ 0,
                      is.na(day_conf)==F && is.na(day_cull_end)==F &&
                        Day_numeric >=day_conf &&
                        Day_numeric <day_cull_end ~ 1,
                      .default=0),
  Sus = case_when(is.na(day_conf)==F && is.na(day_sus)==F &&
                    Day_numeric >=day_sus &&
                    Day_numeric <day_conf ~ 1,
                  is.na(day_conf)==F && is.na(day_cull_end)==F &&
                    Day_numeric >=day_conf &&
                    Day_numeric <day_cull_end ~ 0,
                  .default=0)
) -> cases_daily_long


dat_sf <- st_as_sf(x = cases_daily_long, coords = c("x", "y"), crs = "EPSG:32626")

for (i in 1:max(dat_sf$Day_numeric)){
  tmp = dat_sf %>% filter(Day_numeric==i & Present==1)
  
  p <- ggplot() +
    geom_sf(data = counties, fill = "grey", color = "white",alpha=0.9)+
    geom_sf(data=tmp,aes(color=species),size=1.5)+
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.position="none")
  
  assign(paste0("p", i), p)
} #red = chicken, green = duck

plots <- mget(paste0("p", 1:max(dat_sf$Day_numeric)))
plot_combined = plot_grid(plotlist = plots, ncol = 5,labels=1:max(dat_sf$Day_numeric),align="hvlb")
plot_combined
ggsave("st_plot.png",plot_combined,width=30,height=25,unit="cm",dpi=1000)


