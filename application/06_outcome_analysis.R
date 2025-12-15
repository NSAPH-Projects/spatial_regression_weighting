library(sbw)
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(rnaturalearth)
library(geosphere)
library(gridExtra)
library(Matrix)
library(patchwork)
library(xtable)
library(cowplot)
library(MASS)
library(zoo)
library(duckdb)
library(duckplyr)
library(SuperLearner)
library(geoR)
library(stringr)
library(ggnewscale)
library(data.table)
library(png)
library(scales)
source('funcs.R')

con <- DBI::dbConnect(duckdb())

lab_share <- "" # Fill in with hidden NSAPH lab path on ReD

###################### Step 1: Extract individual level birth claims. ######################

########## 2018 ###########

ip_subset_2018 <- dbGetQuery(con, "
  SELECT  i.clm_id, i.bene_id, i.diagnoses, i.procedures, i.mc_plan_id,
  i.srvc_bgn_dt, i.srvc_end_dt, i.admsn_dt, i.admsn_hr, i.dschrg_dt, i.dschrg_hr,
  b.zcta, b.birth_dt, b.bene_state_cd
  FROM read_parquet(['data/to_data/ip_header/ip_header_2018.parquet']) AS i
  INNER JOIN (
    SELECT bene_id, zcta, birth_dt, bene_state_cd
    FROM read_parquet(['data/to_data/denom/denom_2018.parquet'])
  ) AS b
  ON i.bene_id = b.bene_id
  WHERE i.CLM_TYPE_CD NOT IN ('2', '4', 'B', 'D', 'V', 'X')
  AND EXISTS (
      SELECT 1
      FROM UNNEST(i.diagnoses) AS d(code)
      WHERE REGEXP_MATCHES(code, '^Z38')
  )
")

setDT(ip_subset_2018)              # convert by reference
ip_subset_2018[, clm_id := NULL]   # drop column in-place (like select(-clm_id))
setindex(ip_subset_2018, bene_id)

ip_subset_2018_unique_bene <- ip_subset_2018[
  , .(
    state_cd  = bene_state_cd[1L],
    zcta = zcta[1L],
    srvc_bgn_dt = srvc_bgn_dt[1L],
    srvc_end_dt = srvc_end_dt[1L],
    admsn_dt = admsn_dt[1L],
    admsn_hr = admsn_hr[1L],
    dschrg_dt = dschrg_dt[1L],
    dschrg_hr = dschrg_hr[1L],
    birth_dt = birth_dt[1L],
    diagnoses = list(unlist(diagnoses, recursive = TRUE, use.names = FALSE)),
    procedures = list(unlist(procedures, recursive = TRUE, use.names = FALSE))
  ),
  by = bene_id
]

# Deduplicate also by srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt
ip_subset_2018_unique_bene <- as.data.table(ip_subset_2018_unique_bene)
ip_subset_2018_unique_bene[, .N, by = .(srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt, zcta, birth_dt)][N > 1]
ip_subset_2018_unique_bene <- unique(ip_subset_2018_unique_bene, 
                                     by = c("srvc_bgn_dt", 
                                            "srvc_end_dt", 
                                            "admsn_dt", 
                                            "dschrg_dt", 
                                            "admsn_hr", 
                                            "zcta", 
                                            "birth_dt"))
nrow(ip_subset_2018_unique_bene) # 1,455,152

ip_subset_2018_unique_bene$small_vulnerable <- sapply(ip_subset_2018_unique_bene$diagnoses, function(sublist) {
  if (any(grepl("^P07.0|^P07.1|^P070|^P071|^P07.2|^P072|^P07.3|^P073|^P05|^Z3A0|^Z3A1|^Z3A2|^Z3A30|^Z3A31|^Z3A32|^Z3A33|^Z3A34|^Z3A35|^Z3A36|^O35", sublist))) TRUE else FALSE
})
mean(ip_subset_2018_unique_bene$small_vulnerable) # 14.2%

ip_subset_2018_unique_bene$congenital <- sapply(ip_subset_2018_unique_bene$diagnoses, function(sublist) {
  any(grepl("^Q", sublist))
})
mean(ip_subset_2018_unique_bene$congenital) # 14.6%

state_counts <- ip_subset_2018_unique_bene %>%
  group_by(state_cd) %>%
  summarize(
    n_births = n(),
    n_svn = sum(small_vulnerable, na.rm = T),
    n_congenital = sum(congenital, na.rm = T),
    pct_svn = n_svn/n_births,
    pct_congenital = n_congenital/n_births
  )
mih <- read.csv(paste0(lab_share, 'data/play-doh/cdc_wonder/cdc_wonder_livebirth_counts.csv'))
mih$State.of.Residence.Code <- str_pad(as.character(mih$State.of.Residence.Code), width = 2, pad = '0')
mih_2018 <- subset(mih, Year == 2018)
state_counts <- left_join(state_counts, mih_2018, by = c('state_cd' = 'State.of.Residence.Code'))
state_counts$error <- (state_counts$n_births-state_counts$Births)/state_counts$Births
print(dplyr::select(state_counts, State.of.Residence, n_births, n_svn, pct_svn, n_congenital, pct_congenital, Births, error),
      n = 51)

######## 2017 ###########

ip_subset_2017 <- dbGetQuery(con, "
  SELECT  i.clm_id, i.bene_id, i.diagnoses, i.procedures, i.mc_plan_id,
  i.srvc_bgn_dt, i.srvc_end_dt, i.admsn_dt, i.admsn_hr, i.dschrg_dt, i.dschrg_hr,
  b.zcta, b.birth_dt, b.bene_state_cd
  FROM read_parquet(['data/to_data/ip_header/ip_header_2017.parquet']) AS i
  INNER JOIN (
    SELECT bene_id, zcta, birth_dt, bene_state_cd
    FROM read_parquet(['data/to_data/denom/denom_2017.parquet'])
  ) AS b
  ON i.bene_id = b.bene_id
  WHERE i.CLM_TYPE_CD NOT IN ('2', '4', 'B', 'D', 'V', 'X')
  AND EXISTS (
      SELECT 1
      FROM UNNEST(i.diagnoses) AS d(code)
      WHERE REGEXP_MATCHES(code, '^Z38')
  )
")

setDT(ip_subset_2017)              # convert by reference
ip_subset_2017[, clm_id := NULL]   # drop column in-place (like select(-clm_id))
setindex(ip_subset_2017, bene_id)

ip_subset_2017_unique_bene <- ip_subset_2017[
  , .(
    state_cd  = bene_state_cd[1L],
    zcta = zcta[1L],
    srvc_bgn_dt = srvc_bgn_dt[1L],
    srvc_end_dt = srvc_end_dt[1L],
    admsn_dt = admsn_dt[1L],
    admsn_hr = admsn_hr[1L],
    dschrg_dt = dschrg_dt[1L],
    dschrg_hr = dschrg_hr[1L],
    birth_dt = birth_dt[1L],
    diagnoses = list(unlist(diagnoses, recursive = TRUE, use.names = FALSE)),
    procedures = list(unlist(procedures, recursive = TRUE, use.names = FALSE))
  ),
  by = bene_id
]

# Deduplicate also by srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt
ip_subset_2017_unique_bene <- as.data.table(ip_subset_2017_unique_bene)
ip_subset_2017_unique_bene[, .N, by = .(srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt, zcta, birth_dt)][N > 1]
ip_subset_2017_unique_bene <- unique(ip_subset_2017_unique_bene, 
                                     by = c("srvc_bgn_dt", 
                                            "srvc_end_dt", 
                                            "admsn_dt", 
                                            "dschrg_dt", 
                                            "admsn_hr", 
                                            "zcta", 
                                            "birth_dt"))
nrow(ip_subset_2017_unique_bene) # 1,487,114

ip_subset_2017_unique_bene$small_vulnerable <- sapply(ip_subset_2017_unique_bene$diagnoses, function(sublist) {
  if (any(grepl("^P07.0|^P07.1|^P070|^P071|^P07.2|^P072|^P07.3|^P073|^P05|^Z3A0|^Z3A1|^Z3A2|^Z3A30|^Z3A31|^Z3A32|^Z3A33|^Z3A34|^Z3A35|^Z3A36|^O35", sublist))) TRUE else FALSE
})
mean(ip_subset_2017_unique_bene$small_vulnerable) # 13.9%

ip_subset_2017_unique_bene$congenital <- sapply(ip_subset_2017_unique_bene$diagnoses, function(sublist) {
  any(grepl("^Q", sublist))
})
mean(ip_subset_2017_unique_bene$congenital) # 13.7%

state_counts <- ip_subset_2017_unique_bene %>%
  group_by(state_cd) %>%
  summarize(
    n_births = n(),
    n_svn = sum(small_vulnerable, na.rm = T),
    n_congenital = sum(congenital, na.rm = T),
    pct_svn = n_svn/n_births,
    pct_congenital = n_congenital/n_births
  )
mih_2017 <- subset(mih, Year == 2017)
state_counts <- left_join(state_counts, mih_2017, by = c('state_cd' = 'State.of.Residence.Code'))
state_counts$error <- (state_counts$n_births-state_counts$Births)/state_counts$Births
print(dplyr::select(state_counts, State.of.Residence, n_births, n_svn, pct_svn, n_congenital, pct_congenital, Births, error), 
      n = 51)

####### 2016 ##########
ip_subset_2016 <- dbGetQuery(con, "
  SELECT  i.clm_id, i.bene_id, i.diagnoses, i.procedures, i.mc_plan_id,
  i.srvc_bgn_dt, i.srvc_end_dt, i.admsn_dt, i.admsn_hr, i.dschrg_dt, i.dschrg_hr,
  b.zcta, b.birth_dt, b.bene_state_cd
  FROM read_parquet(['data/to_data/ip_header/ip_header_2016.parquet']) AS i
  INNER JOIN (
    SELECT bene_id, zcta, birth_dt, bene_state_cd
    FROM read_parquet(['data/to_data/denom/denom_2016.parquet'])
  ) AS b
  ON i.bene_id = b.bene_id
  WHERE i.CLM_TYPE_CD NOT IN ('2', '4', 'B', 'D', 'V', 'X')
  AND EXISTS (
      SELECT 1
      FROM UNNEST(i.diagnoses) AS d(code)
      WHERE REGEXP_MATCHES(code, '^Z38')
  )
")

setDT(ip_subset_2016)              # convert by reference
ip_subset_2016[, clm_id := NULL]   # drop column in-place (like select(-clm_id))
setindex(ip_subset_2016, bene_id)

ip_subset_2016_unique_bene <- ip_subset_2016[
  , .(
    state_cd  = bene_state_cd[1L],
    zcta = zcta[1L],
    srvc_bgn_dt = srvc_bgn_dt[1L],
    srvc_end_dt = srvc_end_dt[1L],
    admsn_dt = admsn_dt[1L],
    admsn_hr = admsn_hr[1L],
    dschrg_dt = dschrg_dt[1L],
    dschrg_hr = dschrg_hr[1L],
    birth_dt = birth_dt[1L],
    diagnoses = list(unlist(diagnoses, recursive = TRUE, use.names = FALSE)),
    procedures = list(unlist(procedures, recursive = TRUE, use.names = FALSE))
  ),
  by = bene_id
]

# Deduplicate also by srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt
ip_subset_2016_unique_bene <- as.data.table(ip_subset_2016_unique_bene)
ip_subset_2016_unique_bene[, .N, by = .(srvc_bgn_dt, srvc_end_dt, admsn_dt, dschrg_dt, zcta, birth_dt)][N > 1]
ip_subset_2016_unique_bene <- unique(ip_subset_2016_unique_bene, 
                                     by = c("srvc_bgn_dt", 
                                            "srvc_end_dt", 
                                            "admsn_dt", 
                                            "dschrg_dt", 
                                            "admsn_hr", 
                                            "zcta", 
                                            "birth_dt"))
nrow(ip_subset_2016_unique_bene) # 1419528

ip_subset_2016_unique_bene$small_vulnerable <- sapply(ip_subset_2016_unique_bene$diagnoses, function(sublist) {
  if (any(grepl("^P07.0|^P07.1|^P070|^P071|^P07.2|^P072|^P07.3|^P073|^P05|^Z3A0|^Z3A1|^Z3A2|^Z3A30|^Z3A31|^Z3A32|^Z3A33|^Z3A34|^Z3A35|^Z3A36|^O35", sublist))) TRUE else FALSE
})
mean(ip_subset_2016_unique_bene$small_vulnerable) # 14.2%

ip_subset_2016_unique_bene$congenital <- sapply(ip_subset_2016_unique_bene$diagnoses, function(sublist) {
  any(grepl("^Q", sublist))
})
mean(ip_subset_2016_unique_bene$congenital) # 14.6%

state_counts <- ip_subset_2016_unique_bene %>%
  group_by(state_cd) %>%
  summarize(
    n_births = n(),
    n_svn = sum(small_vulnerable, na.rm = T),
    n_congenital = sum(congenital, na.rm = T),
    pct_svn = n_svn/n_births,
    pct_congenital = n_congenital/n_births
  )
mih_2016 <- subset(mih, Year == 2016)
state_counts <- left_join(state_counts, mih_2016, by = c('state_cd' = 'State.of.Residence.Code'))
state_counts$error <- (state_counts$n_births-state_counts$Births)/state_counts$Births
print(dplyr::select(state_counts, State.of.Residence, n_births, n_svn, pct_svn, n_congenital, pct_congenital, Births, error),
      n = 51)

###################### Step 2: Merge and aggregate to zip code level ################################

ip_subset_2016_2018_unique_bene <- rbind(ip_subset_2016_unique_bene, 
                                         ip_subset_2017_unique_bene,
                                         ip_subset_2018_unique_bene)
nrow(ip_subset_2016_2018_unique_bene) # Total number of Medicaid-delivered births I'm extracting 

state_counts <- ip_subset_2016_2018_unique_bene %>%
  group_by(state_cd) %>%
  summarize(
    n_births = n(),
    n_svn = sum(small_vulnerable, na.rm = T),
    n_congenital = sum(congenital, na.rm = T),
    pct_svn = n_svn/n_births,
    pct_congenital = n_congenital/n_births
  )
state_counts <- subset(state_counts, !(state_cd %in% c(60, 66, 69, 72, 78)))
mih_agg <- mih %>%
  group_by(State.of.Residence.Code) %>%
  summarize(Births = sum(Births),
            State = State.of.Residence[1])
state_counts <- left_join(state_counts, mih_agg, by = c('state_cd' = 'State.of.Residence.Code'))
state_counts$error <- (state_counts$n_births - state_counts$Births)/state_counts$Births
print(dplyr::select(state_counts, State, state_cd, n_births, n_svn, pct_svn, n_congenital, pct_congenital), n = 54)
write.csv(state_counts, 'export/state_counts.csv', row.names = F)
states_error20 <- state_counts$State[abs(state_counts$error) > 0.2]
states_codes_error20 <- state_counts$state_cd[abs(state_counts$error) > 0.2]

# Aggregate the data to the ZCTA level
ip_subset_2016_2018_zcta <- ip_subset_2016_2018_unique_bene %>%
  group_by(zcta) %>%
  summarise(
    state_cd = state_cd[1L],
    total_svn = sum(small_vulnerable, na.rm = TRUE),
    total_congenital = sum(congenital, na.rm = T),
    n_births = n()
  ) %>%
  ungroup()

# Make a map of svn 
zcta_shp <- st_read(paste0(lab_share, "data/lego/geoboundaries/us_geoboundaries__census/us_shapefile__census/zcta_yearly/us_shapefile__census__zcta_yearly__2016/us_shapefile__census__zcta_yearly__2016.shp"))
ip_subset_2016_2018_zcta_shp <- merge(zcta_shp, ip_subset_2016_2018_zcta, by = 'zcta')

# Plot with raw counts by zip code
ggplot() +
  xlim(-125,-65) +
  ylim(25, 50) +
  # Plot polygons, color representing PM2.5 exposure
  geom_sf(data = ip_subset_2016_2018_zcta_shp, aes(fill = 100*total_svn/n_births), color = NA) +
  scale_fill_viridis_c(
    limits = c(0, 20),
    oob = scales::squish,
    breaks = c(0, 5, 10, 20),
    labels = c("0", "5", "10", ">20")
  ) +
  labs(title = "% SVN births in the Medicaid TAP IP, 2016-2018") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),     
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

ggplot() +
  xlim(-125,-65) +
  ylim(25, 50) +
  # Plot polygons, color representing PM2.5 exposure
  geom_sf(data = ip_subset_2016_2018_zcta_shp, aes(fill = 100*total_congenital/n_births), color = NA) +
  scale_fill_viridis_c(
    limits = c(0, 20),
    oob = scales::squish,
    breaks = c(0, 5, 10, 20),
    labels = c("0", "5", "10", ">20")
  ) +
  labs(title = "% Congenital births in the Medicaid TAP IP, 2016-2018") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),     
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

# Plot with censoring to ensure CMS compliance
small_zips <- ip_subset_2016_2018_zcta_shp %>% filter(n_births < 50)
large_zips <- ip_subset_2016_2018_zcta_shp %>% filter(n_births >= 50)

# Find nearest large ZIP for each small one
nearest_idx <- st_nearest_feature(small_zips, large_zips)

# Map small ZIPs to their nearest large ZIP
mapping <- data.frame(
  small_zip = small_zips$zcta,
  nearest_large_zip = large_zips$zcta[nearest_idx]
)
# Merge
zip_summary_merged <- ip_subset_2016_2018_zcta_shp %>%
  filter(state_cd != '01' & state_cd != '08'& state_cd != '44') %>%
  #filter(!(state_cd %in% states_codes_error50)) %>% 
  left_join(mapping, by = c("zcta" = "small_zip")) %>%
  mutate(zip_group = ifelse(is.na(nearest_large_zip), zcta, nearest_large_zip)) %>%
  group_by(zip_group) %>%
  summarize(
    n_births = sum(n_births),
    total_svn = sum(total_svn),
    pct_svn = 100 * total_svn / n_births,
    total_congenital = sum(total_congenital),
    pct_congenital = 100 * total_congenital / n_births
  )

head(zip_summary_merged)
# Plot % SVN by ZCTA, censoring
# Additionally remove 
gsvn <- ggplot() +
  xlim(-125,-65) +
  ylim(25, 50) +
  # Plot polygons, color representing PM2.5 exposure
  geom_sf(data = zip_summary_merged, aes(fill = pct_svn), color = NA) +
  scale_fill_viridis_c(
    name = '% svn',
    limits = c(0, 20),
    oob = scales::squish,
    breaks = c(0, 5, 10, 20),
    labels = c("0", "5", "10", ">20")
  ) +
  labs(title = "% Small vulnerable newborns in the Medicaid TAF Inpatient File, 2016-2018",
       subtitle = "ZCTAs with fewer than 50 live births were merged with the nearest ZCTA \nwith more than 50 live births in accordance with CMS small-cell suppression policy") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_blank(),     
        axis.ticks = element_blank(),
        legend.key.height = unit(30, "points"))


gcong <- ggplot() +
  xlim(-125,-65) +
  ylim(25, 50) +
  # Plot polygons, color representing PM2.5 exposure
  geom_sf(data = zip_summary_merged, aes(fill = pct_congenital), color = NA) +
  scale_fill_viridis_c(
    name = '% congenital',
    limits = c(0, 20),
    oob = scales::squish,
    breaks = c(0, 5, 10, 20),
    labels = c("0", "5", "10", ">20")
  ) +
  labs(title = "% Newborns with congenital anomalies in the Medicaid TAF Inpatient File, 2016-2018",
       subtitle = "ZCTAs with fewer than 50 live births were merged with the nearest ZCTA \nwith more than 50 live births in accordance with CMS small-cell suppression policy") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_blank(),     
        axis.ticks = element_blank(),
        legend.key.height = unit(30, "points"))
png('export/svn.png',
    width = 2000,
    height = 1200,
    res = 220)
plot_with_insets(gsvn)
dev.off()

png('export/cong.png',
    width = 2000,
    height = 1200,
    res = 220)
plot_with_insets(gcong)
dev.off()

################ Step 3: Merge with exposure and confounder data ###################

load('data/preprocessed_superfunds.RData')
n <- nrow(buffers)
buffers <- st_transform(buffers, 5070)

us_outline <- ne_countries(scale = "medium", country = "United States of America", 
                           returnclass = "sf")

ip_subset_2016_2018_zcta_shp$id <- 1:nrow(ip_subset_2016_2018_zcta_shp)
# Transform to equal area crs
ip_subset_2016_2018_zcta_shp <- st_transform(ip_subset_2016_2018_zcta_shp, 5070)

# Compute intersections between buffers and zctas
buffers_union <- st_union(buffers)
ip_subset_2016_2018_zcta_shp <- st_filter(ip_subset_2016_2018_zcta_shp, 
                                          buffers_union, .predicate = st_intersects)
intersections <- st_intersection(buffers, ip_subset_2016_2018_zcta_shp)

# Compute area of each intersection
intersections <- intersections %>% 
  mutate(int_area = as.numeric(st_area(.)))

ip_subset_2016_2018_zcta_shp <- ip_subset_2016_2018_zcta_shp %>% 
  mutate(zcta_area = as.numeric(st_area(.)))

# # Join zcta total area to intersections
intersections <- intersections %>%
  left_join(ip_subset_2016_2018_zcta_shp %>% st_set_geometry(NULL) %>%
              dplyr::select(id, zcta_area), by = "id")

# Calculate the area fraction for each intersection
intersections <- intersections %>% 
  mutate(area_fraction = int_area / zcta_area)

# Estimate births within each intersection, assuming uniform distribution
intersections <- intersections %>% 
  mutate(
    svn_int = total_svn * area_fraction,
    congenital_int = total_congenital * area_fraction,
    num_births = n_births * area_fraction
  )

# Aggregate by buffer (assuming buffer identifier S_EPA_I)
buffers_outcome <- intersections %>%
  group_by(Site_EPA_ID) %>%
  summarise(
    svn_births = sum(svn_int, na.rm = TRUE),
    congenital_births = sum(congenital_int, na.rm = TRUE),
    total_births = sum(num_births, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    svn_percent = 100*svn_births/total_births,
    congenital_percent = 100*congenital_births/total_births,
    total_births = total_births
  )

summary(buffers_outcome$svn_percent)
summary(buffers_outcome$congenital_percent)

buffers_outcome <- left_join(st_drop_geometry(buffers_outcome %>%
                                                dplyr::select(Site_EPA_ID,
                                                              total_births,
                                                              svn_percent,
                                                              congenital_percent)),
                             buffers,
                             by = c('Site_EPA_ID' = 'Site_EPA_ID'))
sum(buffers_outcome$total_births, na.rm = T)
write.csv(sum(buffers_outcome$total_births, na.rm = T), 'export/totalbirths.csv', row.names = F)

################ Step 4: Preliminary exploration for analysis ###################

buffers_outcome <- subset(buffers_outcome, total_births >= 10)
# Drop Alabama, Colorado, Rhode Island because of reporting anomalies
buffers_outcome$state <- substr(buffers_outcome$Site_EPA_ID, 1,2)
buffers_outcome <- subset(buffers_outcome, state != 'AL' & state != 'CO' & state != 'RI')
mod <- lm(svn_percent ~ Z + Site_Score + 
                               population_density + 
                               percent_hispanic + 
                               percent_black + 
                               percent_asian + 
                              percent_indigenous +
                               percent_renter_occupied +
                               metal + voc + water + solid + gas +
                              median_household_income + 
                                median_house_value + 
                               percent_poverty + 
                               percent_high_school_grad + 
                               scale(median_year_structure_built),
          data = buffers_outcome)
summary(mod)$coefficients[2,]

n <- nrow(buffers_outcome) # 1079
lat <- buffers_outcome$Latitude
long <- buffers_outcome$Longitude
X <- buffers_outcome[c('population_density',
                       'percent_hispanic',
                       'percent_black',
                       'percent_indigenous',
                       'percent_asian',
                       'percent_renter_occupied',
                       'median_household_income',
                       'median_house_value',
                       'percent_poverty',
                       'percent_high_school_grad',
                       'median_year_structure_built',
                       'Site_Score',
                       'metal',
                       'voc',
                       'water',
                       'solid',
                       'gas')]

####################### Step 5: print data characteristics #######################

X <- st_drop_geometry(X)
X <- as.matrix(X)
data_characteristics <- data.frame(
  Variable = colnames(X),
  Mean = sapply(1:ncol(X), function(i) mean(X[,i])),
  SD = sapply(1:ncol(X), function(i) sd(X[,i]))
)
# Print the table with xtable for latex without row numbers
print(xtable(data_characteristics, digits = 3), 
      type = "latex", include.rownames = FALSE)
data_characteristics <- rbind.data.frame(c('Binary treatment', mean(buffers_outcome$Z), sd(buffers_outcome$Z)),
                                         data_characteristics)
data_characteristics <- rbind.data.frame(c('% Small Vulnerable', mean(buffers_outcome$svn_percent), sd(buffers_outcome$svn_percent)),
                                         data_characteristics)
data_characteristics <- rbind.data.frame(c('% Congenital', mean(buffers_outcome$congenital_percent), sd(buffers_outcome$congenital_percent)),
                                         data_characteristics)
write.csv(data_characteristics,  file = 'images/summary_stats.csv', row.names = F)

############ Step 6: prepare objects for spatial regression and spatial weighting methods #######################

X <- cbind.data.frame(Intercept = 1, X) # Add intercept

X[,2:ncol(X)] <- scale(X[,2:ncol(X)])
X <- as.matrix(X)

dmat <- distm(cbind(buffers_outcome$Longitude, 
                    buffers_outcome$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

############ Step 7: Run OLS and spatial regression methods to obtain estimates and CIs #######################

outcomes <- cbind.data.frame(Y_svn = buffers_outcome$svn_percent,
                             Y_cong = buffers_outcome$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP')

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')
cis_upper <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_upper) <- methods
rownames(cis_upper) <- c('Small_Vulnerable', 'Congenital')
cis_lower <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_lower) <- methods
rownames(cis_lower) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome$Z,
                      Vre = Vre,
                      Vcar = Vcar,
                      Vgp = Vgp,
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      method = method,
                      lat = lat, 
                      long = long,
                      boot = T)
    
    # Store the estimates
    print(c(fit$est,fit$boot_sd))
    ests[i,j] <- fit$est
    cis_upper[i,j] <- fit$boot_q975
    cis_lower[i,j] <- fit$boot_q025
  }
}

save(ests, cis_lower, cis_upper, file = 'coefficient_estimates/baseline_ests.RData')

############ Step 8: Spatial weighting estimates varying the number of eigenvectors #######################

neigen <- 40
SW_att_ests <- data.frame('neigen' = 1:neigen,
                          'SW_att_svn' = NA,
                          'SW_att_cong' = NA,
                          'ess' = NA)
tols <- vector(mode = "list", length = 3)
names(tols) <- c('RE', 'CAR', 'GP')
tols$RE <- rep(1e-2, neigen)
tols$CAR <- rep(1e-2, neigen)
tols$GP <- rep(1e-2, neigen)

for (i in 1:nrow(SW_att_ests)){
  print(paste('neigen', i, 'of', nrow(SW_att_ests)))
  neigen <- SW_att_ests$neigen[i]
  SW_att_fit_svn <- fit_method(Y = buffers_outcome$svn_percent,
                               X,
                               Z = buffers_outcome$Z,
                               Vre = scale(Vre),
                               Vcar = scale(Vcar),
                               Vgp = scale(Vgp),
                               Sigmainvre = Sigmainvre,
                               Sigmainvcar = Sigmainvcar,
                               Sigmainvgp = Sigmainvgp,
                               tols = tols,
                               method = 'SW',
                               neigen = neigen,
                               boot = F)
  SW_att_ests$SW_att_svn[i] <- SW_att_fit_svn$est
  SW_att_ests$ess[i] <- compute_ess(SW_att_fit_svn$weights)
  
  # Repeat for congenital 
  SW_att_fit_cong <- fit_method(Y = buffers_outcome$congenital_percent,
                                X,
                                Z = buffers_outcome$Z,
                                Vre = scale(Vre),
                                Vcar = scale(Vcar),
                                Vgp = scale(Vgp),
                                Sigmainvre = Sigmainvre,
                                Sigmainvcar = Sigmainvcar,
                                Sigmainvgp = Sigmainvgp,
                                tols = tols,
                                method = 'SW',
                                neigen = neigen,
                                boot = F)
  SW_att_ests$SW_att_cong[i] <- SW_att_fit_cong$est
 
}
plot(SW_att_ests$neigen, SW_att_ests$SW_att_svn, type = 'l')
plot(SW_att_ests$neigen, SW_att_ests$SW_att_cong)

save(SW_att_ests, file = 'coefficient_estimates/SW_att_ests.RData')

# Plot estimates and cis vs neigen
ests_long <- SW_att_ests %>%
  pivot_longer(-neigen, names_to = "outcome", values_to = "att") %>%
  mutate(outcome = str_remove(outcome, "^SW_att_"))

plot_df <- ests_long %>%
  arrange(outcome, neigen) %>%
  dplyr::group_by(outcome) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(SW_att_ests %>% dplyr::select(neigen, ess), by = "neigen")

att_lim <- c(-1, 0.3)
ess_lim <- c(450, 575)

b <- diff(att_lim) / diff(ess_lim)  # (0.2 - (-1)) / (400 - 0) = 1.2/400
a <- att_lim[1] - b * ess_lim[1]    # = -1

png('images/spatial_weighting_estimates.png',
    width = 1500, height = 1500, res = 250)
ggplot(plot_df, aes(x = neigen)) +
  # ATT lines (unified legend: color = outcome, linetype = outcome)
  geom_line(
    aes(y = att,
        color = outcome,
        linetype = outcome),
    linewidth = 1
  ) +
  # ESS line (adds a single ESS entry to the same legend)
  geom_line(
    data = SW_att_ests,
    aes(x = neigen,
        y = a + b * ess,
        color = "ESS",
        linetype = "ESS"),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "svn"  = "steelblue",
      "cong" = "firebrick",
      "ESS"  = "darkgreen"
    ),
    labels = c(
      "svn"  = "Small Vulnerable Newborns",
      "cong" = "Congenital Anomalies",
      "ESS"  = "ESS"
    )
  ) +
  scale_linetype_manual(
    name = "",
    values = c(
      "svn"  = "solid",
      "cong" = "solid",
      "ESS"  = "dashed"
    ),
    labels = c(
      "svn"  = "Small Vulnerable Newborns",
      "cong" = "Congenital Anomalies",
      "ESS"  = "ESS"
    )
  ) +
  scale_y_continuous(
    name = "ATT Estimate",
    limits = att_lim,
    expand = expansion(mult = 0),
    sec.axis = sec_axis(~ (. - a) / b, name = "ESS")
  ) +
  labs(x = "Number of Latent Covariates", title = "") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
dev.off()

############ Step 9: Spatial weighting estimates varying the balancing threshold #######################

k <- 10
threshs <- exp(seq(log(1e-5), log(1e0), length.out = 100))
svn_ests <- rep(NA, length(threshs))
cong_ests <- rep(NA, length(threshs))
esss <- rep(NA, length(threshs))
for (i in 1:length(threshs)){
  thresh <- threshs[i]
  tols$RE[1:k] <- thresh
  tols$CAR[1:k] <- thresh
  tols$GP[1:k] <- thresh
  SW_att_fit_svn <- fit_method(
    Y = buffers_outcome$svn_percent, 
    X,
    Z = buffers_outcome$Z,
    Vre = scale(Vre),
    Vcar = scale(Vcar),
    Vgp = scale(Vgp),
    Sigmainvre = Sigmainvre,
    Sigmainvcar = Sigmainvcar,
    Sigmainvgp = Sigmainvgp,
    tols = tols,
    method = 'SW',
    neigen = k,
    boot = F)
  svn_ests[i] <- SW_att_fit_svn$est
  esss[i] <- compute_ess(SW_att_fit_svn$weights)   
  
  SW_att_fit_cong <- fit_method(
    Y = buffers_outcome$congenital_percent, 
    X,
    Z = buffers_outcome$Z,
    Vre = scale(Vre),
    Vcar = scale(Vcar),
    Vgp = scale(Vgp),
    Sigmainvre = Sigmainvre,
    Sigmainvcar = Sigmainvcar,
    Sigmainvgp = Sigmainvgp,
    tols = tols,
    method = 'SW',
    neigen = k,
    boot = F)
  cong_ests[i] <- SW_att_fit_cong$est
}

plot_df <- tibble(
  threshold = threshs,
  svn = svn_ests,
  cong = cong_ests
) |>
  tidyr::pivot_longer(-threshold, names_to = "outcome", values_to = "att") |>
  dplyr::mutate(outcome = dplyr::recode(
    outcome,
    svn = "Small Vulnerable Newborns",
    cong = "Congenital Anomalies"
  ))
att_lim <- c(-1, 0.3)
ess_lim <- c(525, 575)
b <- diff(att_lim) / diff(ess_lim)        # slope of map
a <- att_lim[1] - b * ess_lim[1]          # intercept'
png('images/balancing_threshold.png',
    width = 1500, height = 1500, res = 250)
ggplot(plot_df, aes(x = threshold)) +
  # ATT lines
  geom_line(aes(y = att, color = outcome, linetype = outcome), linewidth = 1.1) +
  # ESS line (adds a single "ESS" legend entry)
  geom_line(
    data = tibble(threshold = threshs, ess = esss),
    aes(x = threshold, y = a + b * ess, color = "ESS", linetype = "ESS"),
    linewidth = 0.9, inherit.aes = FALSE
  ) +
  scale_x_log10(
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0),
    labels = scales::scientific_format(digits = 1)
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "Small Vulnerable Newborns" = "steelblue",
      "Congenital Anomalies"      = "firebrick",
      "ESS"                       = "darkgreen"
    ),
    breaks = c("ESS", "Congenital Anomalies", "Small Vulnerable Newborns")
  ) +
  scale_linetype_manual(
    name = "",
    values = c(
      "Small Vulnerable Newborns" = "solid",
      "Congenital Anomalies"      = "solid",
      "ESS"                       = "dashed"
    ),
    breaks = c("ESS", "Congenital Anomalies", "Small Vulnerable Newborns")
  ) +
  scale_y_continuous(
    name = "ATT Estimate",
    limits = att_lim,
    expand = expansion(mult = 0),
    sec.axis = sec_axis(~ (. - a) / b, name = "ESS")
  ) +
  labs(x = "Latent Covariate Balancing Threshold", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
dev.off()

################## Step 10: Final spatial weighting estimates and CIs #######################

# Final coefficient estimates
k <- 10 
tols$RE[1:k] <- 1e-2
tols$CAR[1:k] <- 1e-2
tols$GP[1:k] <- 1e-2
SW_att_fit_svn <- fit_method(
  Y = buffers_outcome$svn_percent, 
  X,
  Z = buffers_outcome$Z,
  Vre = scale(Vre),
  Vcar = scale(Vcar),
  Vgp = scale(Vgp),
  Sigmainvre = Sigmainvre,
  Sigmainvcar = Sigmainvcar,
  Sigmainvgp = Sigmainvgp,
  tols = tols,
  method = 'SW',
  neigen = k,
  boot = T) 

# Print as LaTeX table
SW_att_svn <- SW_att_fit_svn$est
SW_lower_svn <- SW_att_fit_svn$boot_q025
SW_upper_svn <- SW_att_fit_svn$boot_q975
print(c(SW_att_svn, SW_lower_svn, SW_upper_svn))

SW_att_fit_cong <- fit_method(
  Y = buffers_outcome$congenital_percent, 
  X,
  Z = buffers_outcome$Z,
  Vre = scale(Vre),
  Vcar = scale(Vcar),
  Vgp = scale(Vgp),
  Sigmainvre = Sigmainvre,
  Sigmainvcar = Sigmainvcar,
  Sigmainvgp = Sigmainvgp,
  tols = tols,
  method = 'SW',
  neigen = k,
  boot = T
)

# Print as LaTeX table
SW_att_cong <- SW_att_fit_cong$est
SW_lower_cong <- SW_att_fit_cong$boot_q025
SW_upper_cong <- SW_att_fit_cong$boot_q975
print(c(SW_att_cong, SW_lower_cong, SW_upper_cong))

# Add to estimates
ests$SW <- NA
ests$SW[1] <- SW_att_svn
ests$SW[2] <- SW_att_cong
cis_lower$SW <- NA
cis_lower$SW[1] <- SW_lower_svn
cis_lower$SW[2] <- SW_lower_cong
cis_upper$SW <- NA
cis_upper$SW[1] <- SW_upper_svn
cis_upper$SW[2] <- SW_upper_cong

### PRINT ESTIMATES
ests_t      <- t(ests)
cis_lower_t <- t(cis_lower)
cis_upper_t <- t(cis_upper)

# Subset only the outcomes you want to include
keep <- c("Small_Vulnerable", "Congenital")

# Format each cell as: est (lower, upper)
formatted <- matrix(nrow = nrow(ests_t), ncol = length(keep))
rownames(formatted) <- rownames(ests_t)
colnames(formatted) <- gsub("^Y_", "", keep)  # cleaner column names

for (i in seq_along(keep)) {
  y <- keep[i]
  formatted[, i] <- sprintf(
    "%.3f (%.3f, %.3f)",
    ests_t[, y],
    cis_lower_t[, y],
    cis_upper_t[, y]
  )
}

# Convert to data frame for xtable
formatted_df <- as.data.frame(formatted)

# save results
write.csv(formatted_df, file = 'images/formatted_df.csv', row.names = F)

################## Step 11: Weighting diagnostics #######################

w_df <- data.frame(
  weight = SW_att_fit_svn$weights,
  Z      = factor(buffers_outcome$Z,
                  levels = c(0, 1),
                  labels = c("Control (Z = 0)", "Treated (Z = 1)"))
)

# 2.  Effective sample size (overall) ------------------------------------
ess <- compute_ess(SW_att_fit_svn$weights)          # a single number
title_txt <- sprintf("ESS = %.1f", ess) # Histogram of spatial weights  (

# 3.  Histogram coloured by treatment group ------------------------------
png('images/spatial_weights_histogram.png',
    width = 1200,
    height = 1500,
    res = 350)
ggplot(w_df, aes(x = weight, fill = Z)) +
  geom_histogram(alpha = 0.65,
                 position = "identity",
                 bins = 30,
                 colour = "grey30") +
  scale_fill_manual(values = c("orange", "dodgerblue")) +
  labs(x = "Spatial weight",
       y = "Count",
       fill = NULL,
       title = title_txt) +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))
dev.off()
img <- readPNG('images/spatial_weights_histogram.png')
grid::grid.raster(img)

colnames(Vre) <- paste0('Vre', 1:ncol(Vre))
colnames(Vcar) <- paste0('Vcar', 1:ncol(Vcar))
colnames(Vgp) <- paste0('Vgp', 1:ncol(Vgp))

# Balance table
bal <- balance_table(
  w = SW_att_fit_svn$weights,
  X = cbind(
    st_drop_geometry(buffers_outcome[
      c("population_density", "percent_hispanic", "percent_black",
        "percent_asian", "percent_indigenous", 
        "percent_renter_occupied", 
        "median_household_income", "median_house_value",
        "percent_poverty", "percent_high_school_grad",
        "median_year_structure_built",
        "Site_Score", "metal", "voc", "solid", "water", "gas")
    ]),
    Vre[, 1:k],
    Vcar[, 1:k],
    Vgp[, 1:k]
  ),
  Z = buffers_outcome$Z
)

bal_round <- bal %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))

# write bal_round to a table
write.csv(bal_round, 'images/bal_round.csv', row.names = F)

# final print
print(xtable(format(bal_round, nsmall = 2), right = TRUE), include.rownames = F)

# NOW: map of weights!  
buffers_outcome$SWiw <- SW_att_fit_svn$weights
buffer_centroids <- st_centroid(buffers[buffers$Site_EPA_ID %in% buffers_outcome$Site_EPA_ID,])
buffers_merged <- merge(buffer_centroids, 
                        dplyr::select(buffers_outcome, Site_EPA_ID, SWiw), 
                        by = 'Site_EPA_ID')

buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)

buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$SWiw),]
buffers_merged_geo <- buffers_merged_geo[order(buffers_merged_geo$Z, decreasing = T),]

const <- unique(buffers_merged_geo$SWiw[buffers_merged_geo$Z == 1])[1]
gSW <- ggplot() +
  geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
  scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
  geom_sf(
    data  = subset(buffers_merged_geo, Z == 0),
    aes(fill = SWiw, shape = factor(Z)),
    size   = 2,
    stroke = .25,
    #alpha = 0.3,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low        = "dodgerblue",
    mid        = "white",
    high       = "orange",
    midpoint   = 0,
    labels     = function(x)
      sprintf("%.3f", x),
    name       = "Control weights",
    guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
  ) +
  new_scale_fill() +
  geom_sf(
    data  = subset(buffers_merged_geo, Z == 1),
    aes(fill = SWiw, shape = factor(Z)),
    size   = 2,
    stroke = .25,
    alpha = 0.5,
    colour = "black"
  ) +
  scale_fill_gradient2(
    low        = "dodgerblue",
    mid = "dodgerblue",
    high = "dodgerblue",
    midpoint   = const,
    limits     = c(const - 1e-4, const + 1e-4),
    breaks     = const,
    labels     = sprintf("%.3f", const),
    name       = "Treated weights",
    guide      = guide_colorbar(order = 3, barheight = unit(0.25, "cm"))  # <<--- shorter
  ) +
  
  labs(shape = "Treatment") + # title = "Spatial Weights", 
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    plot.title   = element_text(hjust = .5),
    legend.position = "right"
  )
png('images/spatial_weights_map.png',
    width = 2000,
    height = 1200,
    res = 220)
plot_with_insets(gSW)
dev.off()
img <- readPNG('images/spatial_weights_map.png')
grid::grid.raster(img)

# PLOT THE Latent COVARIATES
gs <- list()
for (i in 1:k){
  vrei <- Vre[, i]
  buffers_merged <- cbind(buffer_centroids, vrei = vrei)
  buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)
  buffers_merged_geo_reordered <- buffers_merged_geo[order(buffers_merged_geo$vrei, decreasing = T),]
  
  gi <- ggplot() +
    scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
    geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
    geom_sf(
      data  = buffers_merged_geo_reordered,
      aes(fill = vrei, shape = factor(Z)),
      size   = 2,
      stroke = 0,
      colour = "black"
    ) +
    # virids
    scale_fill_viridis_c(
      labels     = function(x)
        sprintf("%.3f", x),
      name       = paste0("Vre", i),
      guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
    ) +
    theme_minimal() +
    # title it with latex \bm{v}_i^{CAR}
    labs(title = bquote(bold(V)[.(i)]^{RE})) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = .5),
      legend.position = "none"
    )
  gs[[i]] <- plot_with_insets(gi)
}

for (i in 1:k){
  vcari <- Vcar[, i]
  buffers_merged <- cbind(buffer_centroids, vcari = vcari)
  buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)
  buffers_merged_geo_reordered <- buffers_merged_geo[order(buffers_merged_geo$vcari),]
  
  gi <- ggplot() +
    scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
    geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
    geom_sf(
      data  = buffers_merged_geo_reordered,
      aes(fill = vcari, shape = factor(Z)),
      size   = 2,
      stroke = 0,
      colour = "black"
    ) +
    # virids
    scale_fill_viridis_c(
      labels     = function(x)
        sprintf("%.3f", x),
      name       = paste0("Vcar", i),
      guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
    ) +
    theme_minimal() +
    # title it with latex \bm{v}_i^{CAR}
    labs(title = bquote(bold(V)[.(i)]^{CAR})) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = .5),
      legend.position = "none"
    )
  gs[[k+i]] <- plot_with_insets(gi)
}

for (i in 1:k){
  vgpi <- Vgp[, i]
  buffers_merged <- cbind(buffer_centroids, vgpi = vgpi)
  buffers_merged_geo <- st_transform(buffers_merged, crs = 4326)
  buffers_merged_geo_reordered <- buffers_merged_geo[order(buffers_merged_geo$vgpi),]
  
  gi <- ggplot() +
    scale_shape_manual(values = c("0" = 21, "1" = 24), guide  = guide_legend(order = 1)) +
    geom_sf(data = us_outline, fill = NA, color = "black", linetype = "solid") +
    geom_sf(
      data  = buffers_merged_geo_reordered,
      aes(fill = vgpi, shape = factor(Z)),
      size   = 2,
      stroke = 0,
      colour = "black"
    ) +
    # virids
    scale_fill_viridis_c(
      labels     = function(x)
        sprintf("%.3f", x),
      name       = paste0("Vgp", i),
      guide      = guide_colorbar(order = 2, barheight = unit(1.0, "cm"))   # <<--- taller
    ) +
    theme_minimal() +
    # title it with latex \bm{v}_i^{gp}
    labs(title = bquote(bold(V)[.(i)]^{GP})) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = .5),
      legend.position = "none"
    )
  gs[[2*k+i]] <- plot_with_insets(gi)
}

png('images/hidden_covariates.png',
    width = 2000, # 2000 # 4000
    height = 900,
    res = 170)
grid.arrange(grobs = gs, nrow = 3)
dev.off()
img <- readPNG('images/hidden_covariates.png')
grid::grid.raster(img)

################## Step 12: Sensitivity analysis varying spatial weighting specification ####################

# 3: Variations
df_sensitivity <- data.frame(model = c('SW_alg', 
                                       'NonlinX_man'),
                             svn = rep(NA, 2),
                             cong = rep(NA, 2))
# Sensitivity 1: SW with algorithm
t_ind <- buffers_outcome$Z
bal_cov <- cbind.data.frame(X,
                            scale(Vre[,1:k]),
                            scale(Vcar[,1:k]), 
                            scale(Vgp[,1:k])) 
colnames(bal_cov) <- paste0('X', 1:ncol(bal_cov))
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal_gri <- c(1e-04,
            0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) # grid of tuning parameters
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_alg <- T # tuning algorithm in Wang and Zubizarreta (2020) used for automatically selecting the degree of approximate covariates balance.
bal$bal_sam <- 1000
bal$bal_std <- 'group'
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_alg_att_svn <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$svn_percent[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$svn_percent[buffers_outcome$Z == 0])
df_sensitivity$svn[1] <- SW_alg_att_svn 
SW_alg_att_congenital <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$congenital_percent[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$congenital_percent[buffers_outcome$Z == 0])
df_sensitivity$cong[1] <- SW_alg_att_congenital 

# Sensitivity 2: Nonlinear only X, manual
t_ind <- buffers_outcome$Z
vars   <- colnames(X[,-1])
sq     <- paste0("I(", vars, "^2)", collapse = " + ")
full_f <- as.formula(paste("~ (.)^2 +", sq))
Xquad <- model.matrix(full_f, data = as.data.frame(X[,-1]))
bal_cov <- cbind.data.frame(Xquad,
                            scale(Vre[,1:k]),
                            scale(Vcar[,1:k]), 
                            scale(Vgp[,1:k])) 
data_frame <- as.data.frame(cbind(t_ind, bal_cov))
t_ind <- "t_ind"
bal <- list()
bal$bal_cov <- colnames(bal_cov)[-1]
bal$bal_std <- 'manual'
bal$bal_alg <- F
bal$bal_tol <- c(rep(0.001, ncol(X)-1), # exact balance on linear terms
                 rep(0.01,ncol(bal_cov)-k*3-ncol(X)), # smallest multiple of 0.01 that would result in consistent constraints
                 tols$RE[1:k], 
                 tols$CAR[1:k], 
                 tols$GP[1:k]) 
sbwatttun_object = sbw(dat = data_frame, ind = t_ind, bal = bal, 
                       sol = list(sol_nam = "quadprog"), 
                       par = list(par_est = "att", par_tar = NULL))
SW_nonlin_att_svn <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$svn_percent[buffers_outcome$Z == 1]) -
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$svn_percent[buffers_outcome$Z == 0])
# Print ATT estimate
df_sensitivity$svn[2] <-SW_nonlin_att_svn
SW_nonlin_att_congenital <- sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 1]*buffers_outcome$congenital_percent[buffers_outcome$Z == 1]) - 
  sum(sbwatttun_object$dat_weights$sbw_weights[buffers_outcome$Z == 0]*buffers_outcome$congenital_percent[buffers_outcome$Z == 0])
df_sensitivity$cong[2] <- SW_nonlin_att_congenital

write.csv(df_sensitivity, file = 'images/df_sensitivity.csv')

################## Step 13: Sensitivity analysis varying the timing of treatment ####################
print(states_error20)
buffers_outcome_error_prone <- subset(buffers_outcome, 
                                      state %in% c('AL', 'AR', 'CA', 'CO', 'HI', 'KS', 'LA', 'MD',
                                                   'MA', 'MN', 'MS', 'NE', 'NV', 'NH', 'NJ', 'RI', 
                                                   'VT', 'WV', 'WY')) # 305
n <- nrow(buffers_outcome_error_prone) # 1392
lat <- buffers_outcome_error_prone$Latitude
long <- buffers_outcome_error_prone$Longitude
X <- buffers_outcome_error_prone[c('population_density',
                       'percent_hispanic',
                       'percent_black',
                       'percent_indigenous',
                       'percent_asian',
                       'percent_renter_occupied',
                       'median_household_income',
                       'median_house_value',
                       'percent_poverty',
                       'percent_high_school_grad',
                       'median_year_structure_built',
                       'Site_Score',
                       'metal',
                       'voc',
                       'water',
                       'solid',
                       'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_error_prone$Longitude, 
                    buffers_outcome_error_prone$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_error_prone$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

outcomes <- cbind.data.frame(Y_svn = buffers_outcome_error_prone$svn_percent,
                             Y_cong = buffers_outcome_error_prone$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')
cis_upper <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_upper) <- methods
rownames(cis_upper) <- c('Small_Vulnerable', 'Congenital')
cis_lower <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_lower) <- methods
rownames(cis_lower) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_error_prone$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10, # constraints inconsistent at 15 and 20
                      boot = F)
    
    # Store the estimates
    ests[i,j]  <- fit$est
  }
}
write.csv(ests, file = 'coefficient_estimates/ests_errorprone.csv')

############################# NON ERROR PRONE STATES #############################
buffers_outcome_nonerror_prone <- subset(buffers_outcome, 
                                         !(Site_EPA_ID %in% buffers_outcome_error_prone$Site_EPA_ID)) # 774

n <- nrow(buffers_outcome_nonerror_prone) 
lat <- buffers_outcome_nonerror_prone$Latitude
long <- buffers_outcome_nonerror_prone$Longitude
X <- buffers_outcome_nonerror_prone[c('population_density',
                                   'percent_hispanic',
                                   'percent_black',
                                   'percent_indigenous',
                                   'percent_asian',
                                   'percent_renter_occupied',
                                   'median_household_income',
                                   'median_house_value',
                                   'percent_poverty',
                                   'percent_high_school_grad',
                                   'median_year_structure_built',
                                   'Site_Score',
                                   'metal',
                                   'voc',
                                   'water',
                                   'solid',
                                   'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_nonerror_prone$Longitude, 
                    buffers_outcome_nonerror_prone$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_nonerror_prone$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

outcomes <- cbind.data.frame(Y_svn = buffers_outcome_nonerror_prone$svn_percent,
                             Y_cong = buffers_outcome_nonerror_prone$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')
cis_upper <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_upper) <- methods
rownames(cis_upper) <- c('Small_Vulnerable', 'Congenital')
cis_lower <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(cis_lower) <- methods
rownames(cis_lower) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_nonerror_prone$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10,
                      boot = F)
    
    # Store the estimates
    ests[i,j] <- fit$est
  }
}
#ests <- rbind.data.frame(ests, cis_lower, cis_upper)
write.csv(ests, file = 'coefficient_estimates/ests_nonerrorprone.csv')

# 2: Timing 
# PERFECT CONTROLS
better_controls <- data$Site_EPA_ID[is.na(data$Construction_Completion_Date) & 
                                      (is.na(data$ra_started) | 
                                         data$ra_started > '2015-12-31')]
buffers_outcome_better_controls <- subset(buffers_outcome, buffers_outcome$Z == 1 | buffers_outcome$Site_EPA_ID %in% better_controls)
n <- nrow(buffers_outcome_better_controls) # 394
lat <- buffers_outcome_better_controls$Latitude
long <- buffers_outcome_better_controls$Longitude
X <- buffers_outcome_better_controls[c('population_density',
                                      'percent_hispanic',
                                      'percent_black',
                                      'percent_indigenous',
                                      'percent_asian',
                                      'percent_renter_occupied',
                                      'median_household_income',
                                      'median_house_value',
                                      'percent_poverty',
                                      'percent_high_school_grad',
                                      'median_year_structure_built',
                                      'Site_Score',
                                      'metal',
                                      'voc',
                                      'water',
                                      'solid',
                                      'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_better_controls$Longitude, 
                    buffers_outcome_better_controls$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_better_controls$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

outcomes <- cbind.data.frame(Y_svn = buffers_outcome_better_controls$svn_percent,
                             Y_cong = buffers_outcome_better_controls$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_better_controls$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10,
                      boot = F)
    
    # Store the estimates
    ests[i,j] <- fit$est
  }
}
write.csv(ests, file = 'coefficient_estimates/ests_bettercontrols.csv')

# LATER CONTROLS
later_controls <- data$Site_EPA_ID[is.na(data$Construction_Completion_Date) & 
                                      (is.na(data$ra_started) | 
                                         data$ra_started >= '2005-01-01')]
buffers_outcome_later_controls <- subset(buffers_outcome, buffers_outcome$Z == 1 | buffers_outcome$Site_EPA_ID %in% later_controls)
n <- nrow(buffers_outcome_later_controls) # 460
lat <- buffers_outcome_later_controls$Latitude
long <- buffers_outcome_later_controls$Longitude
X <- buffers_outcome_later_controls[c('population_density',
                                       'percent_hispanic',
                                       'percent_black',
                                       'percent_indigenous',
                                       'percent_asian',
                                       'percent_renter_occupied',
                                       'median_household_income',
                                       'median_house_value',
                                       'percent_poverty',
                                       'percent_high_school_grad',
                                       'median_year_structure_built',
                                       'Site_Score',
                                       'metal',
                                       'voc',
                                       'water',
                                       'solid',
                                       'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_later_controls$Longitude, 
                    buffers_outcome_later_controls$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_later_controls$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

outcomes <- cbind.data.frame(Y_svn = buffers_outcome_later_controls$svn_percent,
                             Y_cong = buffers_outcome_later_controls$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_later_controls$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10,
                      boot = F)
    
    # Store the estimates
    ests[i,j] <- fit$est
  }
}
write.csv(ests, file = 'coefficient_estimates/ests_latercontrols.csv')

# EARLY TREATED
early_treated <- data$Site_EPA_ID[!(is.na(data$ra_started)) & data$ra_started < '2001-01-01']
buffers_outcome_early_treated <- subset(buffers_outcome, buffers_outcome$Z == 0 | buffers_outcome$Site_EPA_ID %in% early_treated)
n <- nrow(buffers_outcome_early_treated) # 1003
lat <- buffers_outcome_early_treated$Latitude
long <- buffers_outcome_early_treated$Longitude
X <- buffers_outcome_early_treated[c('population_density',
                                       'percent_hispanic',
                                       'percent_black',
                                       'percent_indigenous',
                                       'percent_asian',
                                       'percent_renter_occupied',
                                       'median_household_income',
                                       'median_house_value',
                                       'percent_poverty',
                                       'percent_high_school_grad',
                                       'median_year_structure_built',
                                       'Site_Score',
                                       'metal',
                                       'voc',
                                       'water',
                                       'solid',
                                       'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_early_treated$Longitude, 
                    buffers_outcome_early_treated$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_early_treated$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

tols <- calculate_optimal_tolerances(X = X, 
                                     Z = buffers_outcome_early_treated$Z, 
                                     Sigmainvre = Sigmainvre,
                                     Sigmainvcar = Sigmainvcar,
                                     Sigmainvgp = Sigmainvgp,
                                     Vre = Vre,
                                     Ere = Ere,
                                     Vcar = Vcar,
                                     Ecar = Ecar,
                                     Vgp = Vgp,
                                     Egp = Egp)
outcomes <- cbind.data.frame(Y_svn = buffers_outcome_early_treated$svn_percent,
                             Y_cong = buffers_outcome_early_treated$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_early_treated$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10,
                      boot = F)
    
    # Store the estimates
    ests[i,j] <- fit$est
  }
}
write.csv(ests, file = 'coefficient_estimates/ests_earlytreated.csv')

# LATE TREATED 
late_treated <- data$Site_EPA_ID[!(is.na(data$ra_started)) & 
                                   data$ra_started >= '2001-01-01' & 
                                   data$ra_started <= '2015-12-31']
buffers_outcome_late_treated <- subset(buffers_outcome, buffers_outcome$Z == 0 | buffers_outcome$Site_EPA_ID %in% late_treated)

n <- nrow(buffers_outcome_late_treated) # 894
lat <- buffers_outcome_late_treated$Latitude
long <- buffers_outcome_late_treated$Longitude
X <- buffers_outcome_late_treated[c('population_density',
                                     'percent_hispanic',
                                     'percent_black',
                                     'percent_indigenous',
                                     'percent_asian',
                                     'percent_renter_occupied',
                                     'median_household_income',
                                     'median_house_value',
                                     'percent_poverty',
                                     'percent_high_school_grad',
                                     'median_year_structure_built',
                                     'Site_Score',
                                     'metal',
                                     'voc',
                                     'water',
                                     'solid',
                                     'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome_late_treated$Longitude, 
                    buffers_outcome_late_treated$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

rho2 <- 10
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
S <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(S)
Egp <- E$values
Vgp <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvgp <- solve(Sigma)
# add RE eigens to V
statefactor <- factor(buffers_outcome_late_treated$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
S <- bdiag(Jks)
E <- eigen(S)
Ere <- E$values
Vre <- E$vectors
Sigma <- diag(n) + rho2*S
Sigmainvre <- solve(Sigma)
# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
S <- ginv(as.matrix(L))
E <- eigen(S)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)
Sigma <- diag(n) + rho2*S
Sigmainvcar <- solve(Sigma)

outcomes <- cbind.data.frame(Y_svn = buffers_outcome_late_treated$svn_percent,
                             Y_cong = buffers_outcome_late_treated$congenital_percent)
methods <- c('OLS', 'RE', 'CAR', 'GP', 'SW')  

# Create a dataframe to store estimates and CIs for each method
ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
colnames(ests) <- methods
rownames(ests) <- c('Small_Vulnerable', 'Congenital')

for (i in 1:ncol(outcomes)){
  for (j in 1:length(methods)){
    method <- methods[j]
    print(method)
    Y <- outcomes[,i]
    
    # Fit the model
    fit <- fit_method(Y = Y,
                      X,
                      Z = buffers_outcome_late_treated$Z,
                      Vre = scale(Vre),
                      Vcar = scale(Vcar),
                      Vgp = scale(Vgp),
                      Sigmainvre = Sigmainvre,
                      Sigmainvcar = Sigmainvcar,
                      Sigmainvgp = Sigmainvgp,
                      tols = tols,
                      method = method,
                      lat = lat, 
                      long = long,
                      neigen = 10,
                      boot = F)
    
    # Store the estimates
    print(c(fit$est,fit$boot_sd))
    ests[i,j] <- fit$est
   # <- mean(Y[buffers_outcome_late_treated$Z == 1])/(mean(Y[buffers_outcome_late_treated$Z == 1]) - ATT_est)
  }
}
write.csv(ests, file = 'coefficient_estimates/ests_latetreated.csv')

################# ANOTHER SPATIAL SENSITIVITY ANALYSIS #################
n <- nrow(buffers_outcome) # 1392
lat <- buffers_outcome$Latitude
long <- buffers_outcome$Longitude
X <- buffers_outcome[c('population_density',
                                   'percent_hispanic',
                                   'percent_black',
                                   'percent_indigenous',
                                   'percent_asian',
                                   'percent_renter_occupied',
                                   'median_household_income',
                                   'median_house_value',
                                   'percent_poverty',
                                   'percent_high_school_grad',
                                   'median_year_structure_built',
                                   'Site_Score',
                                   'metal',
                                   'voc',
                                   'water',
                                   'solid',
                                   'gas')]

X <- st_drop_geometry(X)
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
X <- as.matrix(X)
X[,2:ncol(X)] <- scale(X[,2:ncol(X)])

dmat <- distm(cbind(buffers_outcome$Longitude, 
                    buffers_outcome$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
dim(dmat)

# Create an empty sparse matrix for 5-nearest neighbors
nn_matrix <- Matrix(0, n, n, sparse = TRUE)

# Loop through each row to find 5-nearest neighbors
for (i in 1:n) {
  # Get indices of the five smallest distances (excluding the diagonal)
  nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:6]
  
  # Set these indices to 1 in the adjacency matrix
  nn_matrix[i, nearest_indices] <- 1
  nn_matrix[nearest_indices, i] <- 1
}

# Create adjacency matrix
adjacency_matrix <- nn_matrix 

# GP
kappa <- 10 
rangec <- 500
phic <- rangec/(2*sqrt(kappa))
Sgp <- matern(u =  dmat, phi = phic, kappa = kappa)
E <- eigen(Sgp)
Egp <- E$values
Vgp <- E$vectors

# RE
statefactor <- factor(buffers_outcome$state)
Jks <- list()
for (k in unique(statefactor)){
  nk <- sum(statefactor == k)
  Jks[[length(Jks) + 1]] <- matrix(1, nrow = nk, ncol = nk)
}
Sre <- bdiag(Jks)
E <- eigen(Sre)
Ere <- E$values
Vre <- E$vectors

# add CAR eigens to V
L <- diag(rowSums(adjacency_matrix)) - adjacency_matrix
Scar <- ginv(as.matrix(L))
E <- eigen(Scar)
Ecar <- Re(E$values)
Vcar <- Re(E$vectors)

rho2s <- exp(seq(log(0.001), log(20), length.out = 100))  
rho_ests <- list()
for (ix in 1:length(rho2s)){
  rho2 <- rho2s[ix]
  Sigma <- diag(n) + rho2*Sgp
  Sigmainvgp <- solve(Sigma)
  Sigma <- diag(n) + rho2*Sre
  Sigmainvre <- solve(Sigma)
  Sigma <- diag(n) + rho2*Scar
  Sigmainvcar <- solve(Sigma)
  
  tols <- calculate_optimal_tolerances(X = X, 
                                       Z = buffers_outcome$Z, 
                                       Sigmainvre = Sigmainvre,
                                       Sigmainvcar = Sigmainvcar,
                                       Sigmainvgp = Sigmainvgp,
                                       Vre = Vre,
                                       Ere = Ere,
                                       Vcar = Vcar,
                                       Ecar = Ecar,
                                       Vgp = Vgp,
                                       Egp = Egp)
  outcomes <- cbind.data.frame(Y_svn = buffers_outcome$svn_percent,
                               Y_cong = buffers_outcome$congenital_percent)
  methods <- c('OLS', 'RE', 'CAR', 'GP')  
  # Create a dataframe to store estimates and CIs for each method
  ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(methods)))
  colnames(ests) <- methods
  rownames(ests) <- c('Small_Vulnerable', 'Congenital')
  
  for (i in 1:ncol(outcomes)){
    for (j in 1:length(methods)){
      method <- methods[j]
      print(method)
      Y <- outcomes[,i]
      
      # Fit the model
      fit <- fit_method(Y = Y,
                        X,
                        Z = buffers_outcome$Z,
                        Vre = Vre,
                        Vcar = Vcar,
                        Vgp = Vgp,
                        Sigmainvre = Sigmainvre,
                        Sigmainvcar = Sigmainvcar,
                        Sigmainvgp = Sigmainvgp,
                        tols = tols,
                        method = method,
                        lat = lat, 
                        long = long,
                        neigen = 20, 
                        boot = F)
      
      # Store the estimates
      ATT_est <- fit$est
      ests[i,j] <- mean(Y[buffers_outcome$Z == 1])/(mean(Y[buffers_outcome$Z == 1]) - ATT_est)
    }
  }
  rho_ests[[ix]] = ests
}
rho_df <- lapply(seq_along(rho_ests), function(ix) {
  ests <- rho_ests[[ix]]
  as.data.frame(ests) |>
    mutate(Outcome = rownames(ests), rho2 = rho2s[ix]) |>
    pivot_longer(cols = colnames(ests),
                 names_to = "Method", values_to = "Estimate")
}) |>
  bind_rows()

# Rename Congenital and Small_Vulnerable to Congenital Anomalies and Small Vulnerable Newborns
rho_df$Outcome[rho_df$Outcome == 'Small_Vulnerable'] = 'Small Vulnerable Newborns'
rho_df$Outcome[rho_df$Outcome == 'Congenital'] = 'Congenital Anomalies'
png('images/spatial_weighting_estimates_rho2.png',
    width = 2750, height = 1500, res = 250)
ggplot(rho_df, aes(x = rho2, y = Estimate, color = Method)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Outcome, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(
    breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
    labels = label_number(accuracy = 0.01, trim = TRUE)
  ) +
  labs(
    x = expression(rho^2),
    y = "ATT estimates"
  ) +
  ylim(0.94, 1.04) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
dev.off()

# Now for spatal weighting
rho2 <- 10
Sigma <- diag(n) + rho2*Sgp
Sigmainvgp <- solve(Sigma)
Sigma <- diag(n) + rho2*Sre
Sigmainvre <- solve(Sigma)
Sigma <- diag(n) + rho2*Scar
Sigmainvcar <- solve(Sigma)
tols <- calculate_optimal_tolerances(X = X, 
                                     Z = buffers_outcome$Z, 
                                     Sigmainvre = Sigmainvre,
                                     Sigmainvcar = Sigmainvcar,
                                     Sigmainvgp = Sigmainvgp,
                                     Vre = Vre,
                                     Ere = Ere,
                                     Vcar = Vcar,
                                     Ecar = Ecar,
                                     Vgp = Vgp,
                                     Egp = Egp)

multfactor <- exp(seq(log(1e-2), log(1e2), length.out = 100))
rho_ests <- data.frame(matrix(NA, nrow = ncol(outcomes), ncol = length(multfactor)))
rownames(ests) <- c('Small_Vulnerable', 'Congenital')

for (ix in 1:length(multfactor)){
  coeff <- multfactor[ix]
  tols_new <- lapply(tols, FUN = function(x) x * coeff)
  for (i in 1:ncol(outcomes)){
      Y <- outcomes[,i]
      
      # Fit the model
      fit <- fit_method(Y = Y,
                        X,
                        Z = buffers_outcome$Z,
                        Vre = Vre,
                        Vcar = Vcar,
                        Vgp = Vgp,
                        Sigmainvre = Sigmainvre,
                        Sigmainvcar = Sigmainvcar,
                        Sigmainvgp = Sigmainvgp,
                        tols = tols_new,
                        method = 'SW',
                        lat = lat, 
                        long = long,
                        neigen = 20, 
                        boot = F)
      
      # Store the estimates
      ests[i,ix] <- fit$est
    }
}

idx <- seq_len(min(70, length(multfactor), ncol(ests)))

df_plot <- rbind(
  data.frame(
    c_value = multfactor[idx],
    Estimate = as.numeric(ests[1, idx]),
    Outcome = "Small Vulnerable Newborns"
  ),
  data.frame(
    c_value = multfactor[idx],
    Estimate = as.numeric(ests[2, idx]),
    Outcome = "Congenital Anomalies"
  )
)

# Define cube-root transformation
cube_root_trans <- trans_new(
  name = "cube_root",
  transform = function(x) x^(1/3),
  inverse = function(x) x^3
)

png('images/spatial_weighting_estimates_SW_balancingthresh.png',
    width = 1800, height = 1500, res = 250)
ggplot(df_plot, aes(x = c_value, y = Estimate, color = Outcome)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(
    trans = cube_root_trans,
    breaks = pretty_breaks(8),
    labels = label_number(accuracy = 0.01, trim = TRUE)
  ) +
  labs(
    x = expression("Latent covariate balancing threshold"~"= c x "~delta[k]),
    y = "ATT estimate",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
dev.off()
