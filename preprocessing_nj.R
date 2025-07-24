library(tidycensus)
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(sp)
library(spdep)
library(gridExtra)
library(tidyverse)
library(arcpullr)
library(stringr)
library(matrixStats)
library(gstat)
library(geosphere)
library(Matrix)

setwd('/Users/sophie/Documents/implied_weights/superfunds/')

################## 1: Exposure Data ###################

# As of June 6, 2024, there were 1,340 Superfund sites in the National Priorities List in the United States.
# Thirty-nine additional sites have been proposed for entry on the list, 
# and 457 sites have been cleaned up and removed from the list

# Read in csv from arcgis data table in url
# https://epa.maps.arcgis.com/home/item.html?id=c2b7cdff579c41bbba4898400aa38815
url <- 'https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Superfund_National_Priorities_List_(NPL)_Sites_with_Status_Information/FeatureServer/0'
data <- get_spatial_layer(url)

data <- data %>% 
  select(Site_Name, 
         Site_Score, 
         Site_EPA_ID, 
         State, 
         City, 
         County, 
         Status, 
         Latitude, 
         Longitude, 
         Proposed_Date,
         Listing_Date,
         Construction_Completion_Date,
         Deletion_Date) # 1840
# Subset to NJ
data <- data %>% filter(State == 'New Jersey') 
nrow(data) # 153
# Subset those whose Deletion Date is NA or after 2005
# Convert Deletion Date to date
data$Deletion_Date <- as.Date(data$Deletion_Date, format = '%m/%d/%Y')
# Subset to the 1611 superfund sites that were not cleaned up before 2000
data <- data %>% 
  filter(Deletion_Date >= '2001-01-01' | is.na(Deletion_Date))
nrow(data)
# Treated units: superfund sites that have been cleaned up between 2001 and 2010
data$Z <- ifelse(data$Deletion_Date <= '2017-12-31' & data$Deletion_Date >= '2001-01-01', 1, 0)
# Turn NAs in Z to 0 (not yet treated)
data$Z[is.na(data$Z)] <- 0
summary(data$Z) # 6% of superfund sites are treated

data_sf <- st_as_sf(data, coords = c('Longitude', 'Latitude'), crs = 4326)

# Perform the spatial join
# data_sf$cluster <- substr(data_sf$Site_EPA_ID, 1, 2)
# data_sf$cluster <- factor(data_sf$cluster)
# data_sf$cluster <- as.numeric(data_sf$cluster)
# data_sf <- data_sf[order(data_sf$cluster),] # IMPORTANT TO ORDER! 

sf_points_projected <- st_transform(data_sf, crs = 5070)

# Create buffers with a radius of 2 kilometers (2000 meters)
buffers <- st_buffer(sf_points_projected, dist = 2000)

buffers_original_crs <- st_transform(buffers, crs = st_crs(data_sf))
# plot(buffers_original_crs['Status'],xlim = c(-80, -70), ylim = c(38, 46), 
#      main = 'Superfund Sites with 5km Buffers')

setwd('shapefiles/')
st_write(buffers_original_crs, "buffers_2k_nj.shp", append = F)
setwd('../')

#################### 2: Confounder data ###################
#census_api_key("ff4b9c35f317b4718d387706005f65987329e82b", install = T, overwrite = T)

# Read in fips_codes
data(fips_codes)
statefps <- unique(fips_codes$state_code)
statenames <- unique(fips_codes$state_name)
fips_codes <- cbind.data.frame(state_code = statefps, state_name = statenames)
# remove the last 6 rows (not states)
fips_codes <- fips_codes[-c((nrow(fips_codes)-5):nrow(fips_codes)),]

setwd('shapefiles/')
all_shapefiles <- list()
# Unzip and read in each state's tract shapefile in 2000 that has the form: https://www2.census.gov/geo/pvs/tiger2010st/08_Colorado/08/tl_2010_08_tract00.zip
baseURL <- "https://www2.census.gov/geo/pvs/tiger2010st/"
for (state in fips_codes$state_code) {
  state_name <- fips_codes$state_name[fips_codes$state_code == state]
  state_name <- gsub(" ", "_", state_name)
  state_name <- gsub(",", "", state_name)
  state_name <- gsub(" ", "", state_name)
  state_name <- paste0(state, "_", state_name)
  state_name <- paste0(state_name, "/", state, "/")
  url <- paste0(baseURL, state_name, "tl_2010_", state, "_tract00.zip")
  print(url)
  download.file(url, destfile = paste0(state, "_tract00.zip"))
  unzip(paste0(state, "_tract00.zip"))
  # Read in shapefile from the unzipped folder
  state_shapefile <- st_read(paste0('tl_2010_', state, "_tract00.shp"))
  all_shapefiles[[state]] <- state_shapefile
}
# Combine all shapefiles into one
all_shapefiles <- do.call(rbind, all_shapefiles)
setwd('../')

# Define the list of variables for the 2000 Census
vars_2000_sf1 <- c(
  total_population = "P001001",           # Total population
  race_denom = "P009001",                 # Total population for race data
  white_population = "P009002",           # White population"
  black_population = "P009003",           # Black population
  indigenous_population = "P009004",      # American Indian and Alaska Native population
  asian_population = "P009005",           # Asian population
  hispanic_denom = "P004001",             # Total population for Hispanic data
  hispanic_population = "P004002"        # Hispanic population
)

vars_2000_sf3 <- c(
  median_household_income = "P053001",    # Median household income in 1999
  median_house_value = "H085001",         # Median house value
  poverty_population = "P087002",         # Population with income below poverty level in 1999
  poverty_denom = "P087001",              # Total population for poverty rate calculation
  high_school_grad_male = "P037011",      # Male Population with HS degree
  some_college_1_male = "P037012",        # Male Population with <1 year of college
  some_college_male = "P037013",          # Male Population with some college
  college_association_male = "P037014",   # Male Population with Associate's degree
  college_bachelor_male = "P037015",      # Male Population with Bachelor's degree
  college_masters_male = "P037016",       # Male Population with Master's degree
  college_prof_male = "P037017",          # Male Population with Professional degree
  college_doctorate_male = "P037018",     # Male Population with Doctorate degree
  high_school_grad_female = "P037028",    # Female Population with HS degree
  some_college_1_female = "P037029",      # Female Population with <1 year of college
  some_college_female = "P037030",        # Female Population with some college
  college_association_female = "P037031",   # Female Population with Associate's degree
  college_bachelor_female = "P037032",      # Female Population with Bachelor's degree
  college_masters_female = "P037033",       # Female Population with Master's degree
  college_prof_female = "P037034",          # Female Population with Professional degree
  college_doctorate_female = "P037035",     # Female Population with Doctorate degree
  high_school_grad_male_denom = "P037002",     # Total male population aged 25 or older with education info
  high_school_grad_female_denom = "P037019",   # Total female population aged 25 or older with education info
  median_year_built = "H035001"           # Median year housing units were built
)

# Get tract-level data from the 2000 Census
# Initialize empty data frame of 10 columns
tract_data_2000_sf1 <- data.frame(matrix(ncol = 10, nrow = 0))
# Loop through states
for (state in statefps) {
  # Get data for each state
  state_data <- get_decennial(
    geography = "tract",
    year = 2000,
    variables = vars_2000_sf1,
    state = state,
    sumfile = 'sf1',
    output = "wide"
  )
  # Append state data to the empty data frame
  tract_data_2000_sf1 <- rbind(tract_data_2000_sf1, state_data)
}

tract_data_2000_sf3 <- data.frame(matrix(ncol = 10, nrow = 0))
# Loop through states
for (state in statefps) {
  # Get data for each state
  state_data <- get_decennial(
    geography = "tract",
    year = 2000,
    variables = vars_2000_sf3,
    state = state,
    sumfile = 'sf3',
    output = "wide"
  )
  # Append state data to the empty data frame
  tract_data_2000_sf3 <- rbind(tract_data_2000_sf3, state_data)
}

tract_data_2000 <- left_join(tract_data_2000_sf1, 
                             tract_data_2000_sf3, 
                             by = c("GEOID", "NAME"))

# Calculate derived variables (percentages)
tract_data_2000 <- tract_data_2000 %>%
  mutate(
    percent_hispanic = hispanic_population/hispanic_denom,
    percent_black = black_population/race_denom,
    percent_white = white_population/race_denom,
    percent_indigenous = indigenous_population/race_denom,
    percent_asian = asian_population/race_denom,
    percent_poverty = poverty_population / poverty_denom,
    percent_high_school_grad = (
      high_school_grad_male + high_school_grad_female + some_college_1_male + some_college_1_female + some_college_male + some_college_female 
      + college_association_male + college_association_female + college_bachelor_male + college_bachelor_female 
      + college_masters_male + college_masters_female + college_prof_male + college_prof_female 
      + college_doctorate_male + college_doctorate_female
    ) / (high_school_grad_male_denom + high_school_grad_female_denom)
  ) %>%
  select(
    GEOID, NAME, total_population, 
    percent_hispanic, percent_black, percent_white, percent_indigenous,
    percent_asian,  
    median_household_income, median_house_value, 
    percent_poverty, percent_high_school_grad, 
    median_year_built
  )

# View the resulting data
head(tract_data_2000)
summary(tract_data_2000)
# Merge tract_data_2000 with all_shapefiles
tract_data_2000 <- left_join(all_shapefiles, 
                             tract_data_2000, 
                             by = c("CTIDFP00" = "GEOID"))
# change nonsensical values to NA
tract_data_2000 <- tract_data_2000 %>% mutate(
  median_house_value = ifelse(median_house_value == 0, NA, median_house_value),
  median_household_income = ifelse(median_household_income == 0, NA, median_household_income),
  median_year_built = ifelse(median_year_built == 0, NA, median_year_built)
)
save(tract_data_2000, file = 'data/tract_data_2000.RData') # 65443 tracts. 

###################################### 2.5: Spatial Intersection ######################################
load('data/tract_data_2000.RData')
# filter to NJ
tract_data_2000 <- tract_data_2000 %>% filter(STATEFP00 == '34')
tract_data_2000$region <- tract_data_2000$COUNTYFP00

# Read in 2k superfund buffers
buffers <- st_read('shapefiles/buffers_2k_nj.shp')

projected_crs <- 5070

# Reproject both datasets
tract_data_2000 <- st_transform(tract_data_2000, crs = projected_crs)
buffers <- st_transform(buffers, crs = projected_crs)

# Precompute buffer areas and add them as a column
buffers <- buffers %>%
  mutate(buffer_area = as.numeric(st_area(.))) # Convert area to numeric
summary(buffers$buffer_area)

tract_data_2000 <- tract_data_2000 %>%
  mutate(tract_area = as.numeric(st_area(.))) # Calculate and store tract area

# Perform spatial intersection (buffers first to retain buffer info)
intersection <- st_intersection(buffers, tract_data_2000)

# Calculate intersection areas and weighted variables
intersection <- intersection %>%
  mutate(
    intersect_area = as.numeric(st_area(.)), # Intersection area in square meters
    population_weight = intersect_area / tract_area, # Fraction of tract in buffer
    weighted_population = population_weight * total_population # Weighted population
  ) %>%
  mutate(across(
    c(
      'percent_hispanic', 'percent_black', 'percent_white', 'percent_indigenous',
      'percent_asian', 'percent_poverty', 'percent_high_school_grad'
    ),
    ~ .x * intersect_area # Weight percentages by intersection area
  )) # Do NOT modify medians here

buffer_averages <- intersection %>%
  group_by(S_EPA_I) %>% # Group by buffer ID
  summarise(
    total_population = sum(weighted_population, na.rm = TRUE), # Total weighted population
    buffer_area = first(buffer_area), # Precomputed buffer area
    population_density = total_population / buffer_area, # Population density
    across(
      c(
        'percent_hispanic', 'percent_black', 'percent_white', 'percent_indigenous',
        'percent_asian', 'percent_poverty', 'percent_high_school_grad'
      ),
      ~ sum(.x, na.rm = TRUE) / sum(intersect_area, na.rm = TRUE) # Area-weighted averages
    ),
    # Weighted median for median_year_built
    median_year_built = weightedMedian(
      x = median_year_built,
      w = intersect_area,
      na.rm = TRUE
    ),
    # Weighted median for median household income
    median_household_income = weightedMedian(
      x = median_household_income,
      w = intersect_area,
      na.rm = TRUE
    ),
    # Weighted median for median house value
    median_house_value = weightedMedian(
      x = median_house_value,
      w = intersect_area,
      na.rm = TRUE
    ),
    # Region of largest intersection area
    region = names(which.max(table(region)))
  )


# View the resulting buffer-level averages
summary(buffer_averages)

# Which sites got dropped?
buffers$S_EPA_I[which(!buffers$S_EPA_I %in% buffer_averages$S_EPA_I)] # none

# Convert population density (ppl/m^2) to ppl./mi^2
buffer_averages$population_density <- buffer_averages$population_density * 2589988.11
save(buffer_averages, file = 'data/buffer_averages_nj.RData')

#################### 3: Merge exposure and confounder data ####################
buffers <- left_join(buffers, st_drop_geometry(buffer_averages), by = 'S_EPA_I')
buffers <- buffers %>% select(S_EPA_I,
                              Z,
                              region,
                              Latitud,
                              Longitd,
                              population_density,
                              percent_hispanic,
                              percent_black,
                              percent_white,
                              percent_indigenous,
                              percent_asian,
                              median_household_income,
                              median_house_value,
                              percent_poverty,
                              percent_high_school_grad,
                              median_year_built)

# Order by cluster
buffers <- buffers[order(buffers$region),]

# 4: Calculate distance matrix and adjacency matrix
dmat <- distm(cbind(buffers$Longitd, buffers$Latitud), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000
#dmat <- dmat + diag(1e-6, nrow = nrow(dmat))

# Calculate a 2 nearest neighbors matrix using dmat
n <- nrow(dmat)

# Create an empty sparse matrix for 2-nearest neighbors
# nn_matrix <- Matrix(0, n, n, sparse = TRUE)
# 
# # Loop through each row to find 2-nearest neighbors
# for (i in 1:n) {
#   # Get indices of the two smallest distances (excluding the diagonal)
#   nearest_indices <- order(dmat[i, ], decreasing = FALSE)[2:3]
#   
#   # Set these indices to 1 in the adjacency matrix
#   nn_matrix[i, nearest_indices] <- 1
#   nn_matrix[nearest_indices, i] <- 1
# }

# Create adjacency matrix: two points are adjacent if they are within 10k of each other
adjacency_matrix <- dmat < 10 #nn_matrix#

# adjacency_matrix <- st_intersects(buffers_original_crs, sparse = F)
# diag(adjacency_matrix) <- 0

# Save design matrix, distance matrix, and adjacency matrix
X <- buffers[c('population_density',
               'percent_hispanic',
               'percent_black',
               'percent_indigenous',
               'percent_asian',
               'median_household_income',
               'median_house_value',
               'percent_poverty',
               'percent_high_school_grad',
               'median_year_built')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
Z <- buffers$Z
clusters <- buffers$region
# Save preprocessed data. 
save(X, Z, dmat, clusters, adjacency_matrix, buffers, file = 'data/preprocessed_superfunds_nj.RData')

########### Import municipal-level outcome #############
# https://www-doh.nj.gov/doh-shad/query/result/birth/BirthBirthCnty/AvgBirthWt.html
# https://www.nj.gov/health/chs/njshad/?utm_source=chatgpt.com
bw <- readxl::read_xlsx('data/nj_doh_birthweight_municipality.xlsx') 
# Drop first column
bw <- bw[,-1]
# Drop the all NA rows
bw <- bw[complete.cases(bw),]
# If municipality contains the word "unknown", drop the row
bw <- bw[!grepl('unknown', bw$`Mother's Municipality`, ignore.case = TRUE),]
# For now ignore the baseline data, just subset to the years 2018-2022
bw <- bw[bw$`5-Year Groups` == '2018-2022',] # 564 municipalities
# Create a new column for county if the county is in parentheses in the municipality name
bw$County <- ifelse(grepl('\\(.*\\)', bw$`Mother's Municipality`),
                    str_extract(bw$`Mother's Municipality`, '\\(.*\\)'),
                    NA)
# But remove the parentheses from the county
bw$County <- str_remove_all(bw$County, '[()]')
# Remove the county and the extra space from the municipality name
bw$`Mother's Municipality` <- str_remove_all(bw$`Mother's Municipality`, '\\(.*\\)')
bw$`Mother's Municipality` <- str_trim(bw$`Mother's Municipality`)
# Convert Lyndhurst Borough to Township
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Lyndhurst Borough',
                                      'Lyndhurst Township',
                                      bw$`Mother's Municipality`)
# Convert Caldwell Boro Township to Caldwell Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Caldwell Boro Township',
                                      'Caldwell Borough',
                                      bw$`Mother's Municipality`)
# Convert Essex Fells Township to Essex Fells Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Essex Fells Township',
                                      'Essex Fells Borough',
                                      bw$`Mother's Municipality`)
# Convert Glen Ridge Boro Township to Glen Ridge Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Glen Ridge Boro Township',
                                      'Glen Ridge Borough',
                                      bw$`Mother's Municipality`)
# Convert Orange City to City of Orange Township
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Orange City',
                                      'City of Orange Township',
                                      bw$`Mother's Municipality`)
# Convert Princeton Borough to Princeton
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Princeton Borough',
                                      'Princeton',
                                      bw$`Mother's Municipality`)
# Convert Lavalette Borough to Lavallette Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Lavalette Borough',
                                      'Lavallette Borough',
                                      bw$`Mother's Municipality`)
# Convert Peapack and Gladstone Borough to Peapack-Gladstone Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Peapack and Gladstone Borough',
                                      'Peapack-Gladstone Borough',
                                      bw$`Mother's Municipality`)
# Convert Avon-By-The-Sea Borough to Avon-by-the-Sea Borough
bw$`Mother's Municipality` <- ifelse(bw$`Mother's Municipality` == 'Avon-By-The-Sea Borough',
                                      'Avon-by-the-Sea Borough',
                                      bw$`Mother's Municipality`)

# Read in NJ Municipality shapefiles
# https://njogis-newjersey.opendata.arcgis.com/datasets/municipal-boundaries-of-nj-hosted-3424/explore
nj_mun <- st_read('shapefiles/NJ_Municipal_Boundaries/NJ_Municipal_Boundaries_3424.shp')
mean(bw$`Mother's Municipality` %in% nj_mun$MUN_LABEL) # 0.99.8
nj_mun$bw <- NA
for (i in 1:nrow(bw)){
  if (bw$`Mother's Municipality`[i] %in% nj_mun$MUN_LABEL){
    # If one match for municipality name, 
    if (is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$`Mother's Municipality`[i])
      nj_mun$bw[j] <- bw$`Average Birth Weight of Infants (grams)`[i] 
    }
    # If two matches for municipality name, choose the one with the same county
    if (!is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$`Mother's Municipality`[i] & nj_mun$COUNTY == toupper(bw$County[i]))
      nj_mun$bw[j] <- bw$`Average Birth Weight of Infants (grams)`[i] 
    }
  }
}
nj_mun$bw <- as.numeric(gsub(',', '', nj_mun$bw))
# Convert bw from g to lbs

png('images/municipality_bw_nj.png', width = 1000, height = 1000, res = 150)
ggplot() +
  geom_sf(data = nj_mun, aes(fill = bw), color = NA) +
  labs(fill = 'Average Birthweight (lbs)') +
  ggtitle('5-Year Average Birthweight by Municipality, 2018-2022') +
  scale_fill_viridis_c() +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

load('data/preprocessed_superfunds_nj.RData')
# For each buffer, assign to it a Y value that is average percent low bw
nj_mun <- st_transform(nj_mun, st_crs(buffers))
intersections <- st_intersection(buffers, nj_mun) # 389 intersections but only 129 EPA IDs... 
# Add area of intersection
intersections$area <- st_area(intersections)
intersections$area <- as.numeric(intersections$area)
buffers_outcome <- intersections %>% 
  group_by(S_EPA_I) %>% # compute weighted mean
  summarise(Y = weighted.mean(bw, w = area, na.rm = T))
summary(buffers_outcome$Y)

# Merge with buffers
buffers_outcome <- left_join(buffers, 
                             st_drop_geometry(buffers_outcome %>% 
                                                select(S_EPA_I,Y)), 
                             by = c('S_EPA_I' = 'S_EPA_I'))
summary(buffers_outcome$Y)
# For some reason a superfund site in Howell NJ was not included in the merge... manually add
bwHowell <- bw$`Average Birth Weight of Infants (grams)`[bw$`Mother's Municipality` == 'Howell Township']
# Remove comma and convert to numeric like above
bwHowell <- as.numeric(gsub(',', '', bwHowell))
buffers_outcome$Y[is.na(buffers_outcome$Y)] <- bwHowell
summary(buffers_outcome$Y)

class(buffers_outcome)
# Save 
save(X, Z, dmat, clusters, adjacency_matrix, buffers_outcome, 
     file = 'data/preprocessed_superfunds_outcome_nj.RData')

model_bw <- lm(Y ~ Z + log(population_density) + 
              percent_hispanic +
              percent_black +
              percent_indigenous +
              percent_asian +
              log(median_household_income) + 
              log(median_house_value) + 
              percent_poverty + 
              percent_high_school_grad,
            data = buffers_outcome)
summary(model_bw)
model_bw$coefficients['Z']
# Convert grams to pounds
model_bw$coefficients['Z'] * 0.00220462

# Add pm25 from buffer_centroids_ordered to buffers_outcome, merging on S_EPA_I
buffers_outcome <- left_join(buffers_outcome, 
                             st_drop_geometry(buffer_centroids_ordered) %>% 
                               select(S_EPA_I, pm25), 
                             by = c('S_EPA_I' = 'S_EPA_I'))
model_bw_pm25 <- lm(Y ~ Z + log(population_density) + 
                 percent_hispanic +
                 percent_black +
                 percent_indigenous +
                 percent_asian +
                 log(median_household_income) + 
                 log(median_house_value) + 
                 percent_poverty + 
                 percent_high_school_grad + 
                   pm25,
               data = buffers_outcome)
# Extract pm25 coefficient
model_bw_pm25$coefficients['pm25']* 0.00220462 # -0.0008857245 
model_bw_pm25$coefficients['Z'] * 0.00220462
buffers <- left_join(buffers, 
                     st_drop_geometry(buffers_outcome) %>% 
                       select(S_EPA_I, Y), 
                     by = c('S_EPA_I' = 'S_EPA_I'))

# sum(l*buffers$Y)
#[1] 2.392094

# Add Site score from buffers$Sit_Scr to buffers_outcome, merging on S_EPA_I
buffers_outcome <- left_join(buffers_outcome, 
                             st_drop_geometry(buffers) %>% 
                               select(S_EPA_I, Sit_Scr), 
                             by = c('S_EPA_I' = 'S_EPA_I'))
model_bw_Sit_Scr <- lm(Y ~ Z + log(population_density) +
                              percent_hispanic +
                              percent_black +
                              percent_indigenous +
                              percent_asian +
                              log(median_household_income) + 
                              log(median_house_value) + 
                              percent_poverty + 
                              percent_high_school_grad + 
                              Sit_Scr,
                            data = buffers_outcome)
# Extract Site Score coefficient
model_bw_Sit_Scr$coefficients['Sit_Scr']* 0.00220462 # -0.000718599
model_bw_Sit_Scr$coefficients['Z'] * 0.00220462

# Compute Moran's I of residuals
residuals <- residuals(model_bw)
# Inverse Distance based weights
coords <- cbind(buffers_outcome$Latitud, buffers_outcome$Longitd)
dmat <- distm(cbind(buffers_outcome$Longitd, buffers_outcome$Latitud), fun = distHaversine)
# Avoid division by zero by adding a small constant to zero distances
distance_weights <- 1/dmat
diag(distance_weights) <- 0
distance_weights_listw <- mat2listw(distance_weights, style = 'W', zero.policy = T)
# Moran test
moran.test(residuals, distance_weights_listw) # suggests weak but significant positive spatial autocorrelation

# 1. Compute the Haversine distance matrix
# Assuming buffers_outcome contains longitude ('Longitd') and latitude ('Latitud')
dmat <- distm(cbind(buffers_outcome$Longitd, buffers_outcome$Latitud), fun = distHaversine)

# 2. Define bins up to 250,000 meters (250 km) with a 25,000-meter (25 km) bin width
bin_width <- 7500  
distance_bins <- seq(0, 150000, by = bin_width)

# 3. Compute semivariance for each bin
semivariance <- sapply(1:(length(distance_bins) - 1), function(k) {
  # Get pairs in the current bin
  lower_bound <- distance_bins[k]
  upper_bound <- distance_bins[k + 1]
  in_bin <- (dmat >= lower_bound & dmat < upper_bound)
  
  # Extract residual pairs
  residual_pairs <- outer(residuals, residuals, `-`)[in_bin]
  
  # Compute semivariance
  if (length(residual_pairs) > 0) {
    return(mean(residual_pairs^2) / 2)
  } else {
    return(NA)  # No pairs in this bin
  }
})

# 4. Calculate midpoints of bins for plotting
bin_midpoints <- (distance_bins[-1] + distance_bins[-length(distance_bins)]) / 2

# 5. Plot the semivariogram
# plot(bin_midpoints / 1000, semivariance, pch = 16, 
#      xlab = "Distance (km)", ylab = "Semivariance", 
#      main = "Semivariogram of Residuals (0–500 km)", type = 'l')


# Create a data frame for ggplot2
semivariogram_df <- data.frame(
  Distance_km = bin_midpoints / 1000,  # Convert meters to kilometers
  Semivariance = semivariance
)

png('images/semivariogram_nj.png', width = 1000, height = 600, res = 150)
# Plot using ggplot2
ggplot(semivariogram_df, aes(x = Distance_km, y = Semivariance)) +
  geom_point(size = 3) +  # Points for each semivariance
  geom_line(linewidth = 1) +   # Line connecting the points
  labs(
    title = "Semivariogram of Residuals (0–150 km)",
    x = "Distance (km)",
    y = "Semivariance"
  ) +
  theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()
