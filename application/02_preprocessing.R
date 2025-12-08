library(dplyr)
library(sf)
library(arcpullr)
library(matrixStats)
library(geosphere)
library(Matrix)
library(mice)
# For web scraping EPA Superfund schedule pages
library(rvest)
library(dplyr)
library(stringr)
library(xml2)
library(parallel)


setwd('../') # Work in parent directory

################## 1: Exposure Data ###################

# Read in csv from arcgis data table in url: Superfund Sites Where you Live EPA ArcGIS
# url <- 'https://epa.maps.arcgis.com/home/item.html?id=c2b7cdff579c41bbba4898400aa38815'
# https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=33cebcdfdd1b4c3a8b51d416956c41f1
url <- 'https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Superfund_National_Priorities_List_(NPL)_Sites_with_Status_Information/FeatureServer/0'
data <- get_spatial_layer(url)

data <- data %>% 
  dplyr::select(Site_Name, 
                SEMS_ID,
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
# What is the average time between listing date and deletion date
data$Listing_Date <- as.Date(data$Listing_Date, format = '%m/%d/%Y')
data$Deletion_Date <- as.Date(data$Deletion_Date, format = '%m/%d/%Y')
data$Construction_Completion_Date <- as.Date(data$Construction_Completion_Date, 
                                             format = '%m/%d/%Y')
data$Proposed_Date <- as.Date(data$Proposed_Date, format = '%m/%d/%Y')

# Subset to those actually listed
data <- data %>% 
  filter(!(is.na(Listing_Date)))
# # Subset to those not deleted before 2001
data <- data %>%
  filter(is.na(Deletion_Date) | Deletion_Date >= '1991-01-01')
nrow(data) 

# Web Scrape to get Remediation Action Started Date
data$link <- paste0('https://cumulis.epa.gov/supercpad/SiteProfiles/index.cfm?fuseaction=second.schedule&id=', 
                    str_pad(data$SEMS_ID, pad = '0', width = 7))
data$ra_started <- NA

start_time <- Sys.time()
data$ra_started <- unlist(
  mclapply(
    data$link,
    safe_get_ra_started,
    mc.cores = 6,   
    mc.preschedule = FALSE
  )
)
end_time <- Sys.time()
end_time - start_time

save(data, file = '../data/superfund_sites_with_ra_started.RData')

# Deal with the 11 NAs manually 
data$ra_started[125] <- 'Not Yet Achieved'
data$ra_started[297] <- '06/30/1992' # There are several operable units, each with separate timelines. taking the earliest remedial action date.
data$ra_started[396] <- 'Not Yet Achieved'
data$ra_started[675] <- 'Not Yet Achieved'
data$ra_started[676] <- 'Not Yet Achieved'
data$ra_started[750] <- 'Not Yet Achieved'
data$ra_started[1025] <- 'Not Yet Achieved'
data$ra_started[1348] <- 'Not Yet Achieved'
data$ra_started[1657] <- 'Not Yet Achieved'
data$ra_started[1670] <- '01/01/1990' # see cleanup timeline https://www.epa.gov/wyckoff-eagle-harbor/site-history
data$ra_started[1706] <- 'Not Yet Achieved'

# Convert 'Not Yet Achieved' and anything that starts with 'Estimated' in ra_started to NA
data$ra_started <- ifelse(grepl('Not Yet Achieved', data$ra_started) | grepl('Estimated', data$ra_started), 
                          NA, data$ra_started)

data$ra_started <- as.Date(data$ra_started, format = '%m/%d/%Y')

sum(is.na(data$ra_started) | data$ra_started >= '1991-01-01') # 1473 had remedial action start after 1991
# 1458
sum(data$Deletion_Date < '2015-12-31' & (is.na(data$ra_started) | data$ra_started >= '1991-01-01'), na.rm = T) # 278 were cleaned up by 2015 and had remedial action start after 1991
# 263
sum((data$Deletion_Date >= '2016-01-01' | is.na(data$Deletion_Date)) & 
      (is.na(data$ra_started) | data$ra_started >= '1991-01-01'), na.rm = T) # 1195 were not yet cleaned up by 2016 and had remedial action start after 1991
# 1195

# Subset to the Superfund sites that had remedial action start after 1991 (or not started yet)
data <- data %>% 
  filter(is.na(ra_started) | ra_started >= '1991-01-01') # remedial action started after 1991
nrow(data)

# data$time_to_deletion <- as.numeric(difftime(data$Deletion_Date, data$Listing_Date, units = "days"))/365.25
# summary(data$time_to_deletion) # 1-41 years, median 13, mean 16 years

########### READ IN CHEMICALS AND MEDIUMS #############
contam <- read.csv('data/site_chemical_medium_flags.csv')
# Not all superfund sites have contaminants reported there
# Merge contam with data
data <- left_join(data, contam, by = c('Site_EPA_ID' = 'S_EPA_I'))

# Treated units: superfund sites that have been cleaned up between 2001 and 2015
data$Z <- ifelse(data$Deletion_Date < '2015-12-31', 1, 0) 
# Turn NAs in Z to 0 (not yet treated)
data$Z[is.na(data$Z)] <- 0
sum(data$Z) # 263
summary(data$Z) # 18% of superfund sites are treated
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326, remove = F)

# https://www.arcgis.com/home/item.html?id=d6e1591d9a424f1fa6d95a02095a06d7
# Site boundaries are polygons representing the footprint of a whole site, 
# defined for purposes of this effort as the sum of all of the 
# Operable Units and the current understanding of the full extent of contamination. 
url <- 'https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/FAC_Superfund_Site_Boundaries_EPA_Public/FeatureServer/0'
boundaries <- get_spatial_layer(url)
boundaries <- st_transform(boundaries, st_crs(data_sf))

# Assign boundaries as geometry of data
idx <- match(data_sf$Site_EPA_ID, boundaries$EPA_ID)   # index of matching boundary per row
has <- !is.na(idx) # 99 percent

g_data <- st_geometry(data_sf)
g_bnd  <- st_geometry(boundaries)

g_data[has] <- g_bnd[idx[has]]
st_geometry(data_sf) <- g_data

# Perform the spatial join
data_sf$cluster <- substr(data_sf$Site_EPA_ID, 1, 2)
data_sf$cluster <- factor(data_sf$cluster)
data_sf$cluster <- as.numeric(data_sf$cluster)
data_sf <- data_sf[order(data_sf$cluster),]  

sf_points_projected <- st_transform(data_sf, crs = 5070)

# Create buffers with a radius of 2 kilometers (2000 meters)
buffers <- st_buffer(sf_points_projected, dist = 2000)

buffers_original_crs <- st_transform(buffers, crs = st_crs(data_sf))

st_write(buffers_original_crs, "shapefiles/buffers_2k.shp", append = F)

#################### 2: 1990 Confounder data ###################
# https://catalog.data.gov/dataset/census-tracts-in-1990/resource/18c7c256-0b30-41e0-887e-5b0b2ce2dc20
tract_1990_shp <- st_read('data/nhgis0001_shape/nhgis0001_shapefile_tl2000_us_tract_1990/US_tract_1990.shp')
# population density, percentage of individuals identifying as Hispanic, Black, Indigenous, or Asian, 
# median household income, median house value, percentage below the poverty line, tenure status,
# percentage with a high school diploma, and the median year of housing construction

# Sourced from IPUMS NHGIS
# Summary File 1
sf1 <- read.csv('data/nhgis0002_csv/nhgis0002_ds120_1990_tract.csv')
sf1 <- dplyr::select(sf1, GISJOIN, 
                     ET1001, 
                     EUY001, EUY002, EUY003, EUY004, EUY005,
                     EU0001,
                     ES1001, ES1002,
                     EST001)
colnames(sf1)[2:ncol(sf1)] <- c(
  'total_population',
  'white_population', 'black_population', 'indigenous_population', 'asian_population', 'other_race',
  'hispanic_population',
  'owner_occupied', 'renter_occupied',
  'median_house_value')
summary(sf1)
 
# Summmary File 3
sf3 <- read.csv('data/nhgis0002_csv/nhgis0002_ds123_1990_tract.csv')
sf3 <- dplyr::select(sf3, GISJOIN,
                     E33001, E33002, E33003, E33004, E33005, E33006, E33007,
                     E4U001,
                     E07001, E07002, E07003, E07004, E07005, E07006, E07007, E07008,
                     E07009, E07010, E07011, E07012,
                     E07013, E07014, E07015, E07016, E07017, E07018, E07019, E07020,
                     E07021, E07022, E07023, E07024,
                     EX8001)


colnames(sf3)[2:ncol(sf3)] <- c(
  'less_than_9th_grade', 'some_HS_no_diploma', 'high_school_grad', 'some_college_no_degree', 'associate_degree', 'bachelor_degree', 'graduate_professional_degree',
  'median_household_income',
  'income_above_poverty_under_5', 'income_above_poverty_5', 'income_above_poverty_6_to_11',
  'income_above_poverty_12_to_17', 'income_above_poverty_18_to_24', 'income_above_poverty_25_to_34',
  'income_above_poverty_35_to_44', 'income_above_poverty_45_to_54',
  'income_above_poverty_55_to_59', 'income_above_poverty_60_to_64', 'income_above_poverty_65_to_74',
  'income_above_poverty_75_and_over',
  'income_below_poverty_under_5', 'income_below_poverty_5', 'income_below_poverty_6_to_11',
  'income_below_poverty_12_to_17', 'income_below_poverty_18_to_24', 'income_below_poverty_25_to_34',
  'income_below_poverty_35_to_44', 'income_below_poverty_45_to_54',
  'income_below_poverty_55_to_59', 'income_below_poverty_60_to_64', 'income_below_poverty_65_to_74',
  'income_below_poverty_75_and_over',
  'median_year_structure_built'
)

# Replace the 0 values in median household income, median house value, and median year structure built with NA
sf1$median_house_value[sf1$median_house_value == 0] <- NA
sf3$median_household_income[sf3$median_household_income == 0] <- NA
sf3$median_year_structure_built[sf3$median_year_structure_built == 0] <- NA

# Merge sf1 and sf3
tract_data_1990 <- left_join(sf1, sf3, by = 'GISJOIN')

# Calculate the derived variables
tract_data_1990 <- tract_data_1990 %>%
  mutate(
    total_race_population = white_population + black_population + indigenous_population + asian_population + other_race,
    percent_hispanic = hispanic_population / total_population,
    percent_black = black_population / total_race_population,
    percent_white = white_population / total_race_population,
    percent_indigenous = indigenous_population / total_race_population,
    percent_asian = asian_population / total_race_population,
    percent_renter_occupied = renter_occupied / (owner_occupied + renter_occupied),
    percent_poverty = (income_below_poverty_under_5 + income_below_poverty_5 + income_below_poverty_6_to_11 +
                         income_below_poverty_12_to_17 + income_below_poverty_18_to_24 + income_below_poverty_25_to_34 +
                         income_below_poverty_35_to_44 + income_below_poverty_45_to_54 +
                         income_below_poverty_55_to_59 + income_below_poverty_60_to_64 + income_below_poverty_65_to_74 +
                         income_below_poverty_75_and_over) /
      (income_above_poverty_under_5 + income_above_poverty_5 + income_above_poverty_6_to_11 +
         income_above_poverty_12_to_17 + income_above_poverty_18_to_24 + income_above_poverty_25_to_34 +
         income_above_poverty_35_to_44 + income_above_poverty_45_to_54 +
         income_above_poverty_55_to_59 + income_above_poverty_60_to_64 + income_above_poverty_65_to_74 +
         income_above_poverty_75_and_over +
         income_below_poverty_under_5 + income_below_poverty_5 + income_below_poverty_6_to_11 +
         income_below_poverty_12_to_17 + income_below_poverty_18_to_24 + income_below_poverty_25_to_34 +
         income_below_poverty_35_to_44 + income_below_poverty_45_to_54 +
         income_below_poverty_55_to_59 + income_below_poverty_60_to_64 + income_below_poverty_65_to_74 +
         income_below_poverty_75_and_over),
    percent_high_school_grad = (high_school_grad + some_college_no_degree + associate_degree + bachelor_degree + graduate_professional_degree) /
      (less_than_9th_grade + some_HS_no_diploma + high_school_grad + some_college_no_degree + associate_degree + bachelor_degree + graduate_professional_degree)) %>%
    dplyr::select(
      GISJOIN, total_population, 
      percent_hispanic, percent_black, percent_white, percent_indigenous,
      percent_asian, percent_renter_occupied,
      median_household_income, median_house_value, 
      percent_poverty, percent_high_school_grad, 
      median_year_structure_built
    )

# View the resulting data
head(tract_data_1990)
summary(tract_data_1990)
# Merge tract_data_1990 with shapefile
tract_data_1990 <- left_join(tract_1990_shp, 
                             tract_data_1990, 
                             by = c("GISJOIN"))
save(tract_data_1990, file = 'data/tract_data_1990.RData') # 60947 tracts 

#################################### Merge with buffers  ####################################

# Read in 2k superfund buffers
buffers <- st_read('shapefiles/buffers_2k.shp')

projected_crs <- 5070

# Reproject both datasets
tract_data_1990 <- st_transform(tract_data_1990, crs = projected_crs)
buffers <- st_transform(buffers, crs = projected_crs)

# Precompute buffer areas and add them as a column
buffers <- buffers %>%
  mutate(buffer_area = as.numeric(st_area(.))) # Convert area to numeric
# Summarize buffer areas in square mile
quantile(buffers$buffer_area/2589988.11, probs = c(0.1,0.9)) # 5 - 23 square miles typically. but biggest is 1000

tract_data_1990 <- tract_data_1990 %>%
  mutate(tract_area = as.numeric(st_area(.))) # Calculate and store tract area

# Perform spatial intersection (buffers first to retain buffer info)
intersection <- st_intersection(buffers, tract_data_1990)

# Calculate intersection areas and weighted variables
intersection <- intersection %>%
  mutate(
    intersect_area = as.numeric(st_area(.)), # Intersection area in square meters
    population_weight = intersect_area / tract_area, # Fraction of tract in buffer
    weighted_population = population_weight * total_population # Estimated pop if population was uniformly distributed
  ) %>%
  mutate(across(
    c(
      'percent_hispanic', 'percent_black', 'percent_white', 'percent_indigenous',
      'percent_asian', 'percent_poverty', 'percent_high_school_grad', 'percent_renter_occupied'
    ),
    ~ .x * intersect_area # Weight percentages by intersection area
  )) # Do NOT modify medians here

buffer_averages <- intersection %>%
  group_by(Site_EPA_ID) %>% # Group by buffer ID
  summarise(
    total_population = sum(weighted_population, na.rm = TRUE), # Total estimated population if uniform
    buffer_area = first(buffer_area), # Precomputed buffer area
    population_density = total_population / buffer_area, # Population density
    across(
      c(
        'percent_hispanic', 'percent_black', 'percent_white', 'percent_indigenous',
        'percent_asian', 'percent_poverty', 'percent_high_school_grad', 'percent_renter_occupied'
      ),
      ~ sum(.x, na.rm = TRUE) / sum(intersect_area, na.rm = TRUE) # Area-weighted averages (previously weighted above)
    ),
    # Weighted median for median_year_structure_built
    median_year_structure_built = weightedMedian(
      x = median_year_structure_built,
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
    )
  )


# View the resulting buffer-level averages
summary(buffer_averages)

# Now n = 1584. Which sites got dropped?
buffers$Site_EPA_ID[which(!buffers$Site_EPA_ID %in% buffer_averages$Site_EPA_ID)] # PR, Guam, Virgin Islands

# Convert population density (ppl/m^2) to ppl./mi^2
buffer_averages$population_density <- buffer_averages$population_density * 2589988.11
save(buffer_averages, file = 'data/buffer_averages.RData')

#################### 3: Merge exposure and confounder data ####################
# rename the column Site_EPA_ID in buffer_averages to S_EPA_I
# buffer_averages <- buffer_averages %>%
#   rename(S_EPA_I = Site_EPA_ID)
buffers <- buffers %>%
  left_join(
    st_drop_geometry(buffer_averages),
    by = c("Site_EPA_ID")
  )
buffers <- buffers %>% dplyr::select(Site_EPA_ID,
                                     Site_Score,
                                     metal,
                                     voc, # nor including pop or pah since always present
                                     gas,
                                     solid,
                                     water,
                                      Z,
                                      cluster,
                                      Latitude,
                                      Longitude,
                                      population_density,
                                      percent_hispanic,
                                      percent_black,
                                      percent_white,
                                      percent_indigenous,
                                      percent_asian,
                                      percent_renter_occupied,
                                      median_household_income,
                                      median_house_value,
                                      percent_poverty,
                                      percent_high_school_grad,
                                      median_year_structure_built)
# Drop the rows with NA census confounder values (PR, Guam, Virgin Islands)
buffers <- buffers[!is.na(buffers$percent_hispanic),] # 1458 --> 1429
# That line wasn't specific to missingness in Hispanic but census confounders more generally
# Order by cluster
buffers <- buffers[order(buffers$cluster),]
# Convert metal, voc, gas, solid, water to factors
buffers$metal <- as.factor(buffers$metal)
buffers$voc <- as.factor(buffers$voc)
buffers$gas <- as.factor(buffers$gas)
buffers$solid <- as.factor(buffers$solid)
buffers$water <- as.factor(buffers$water)

# Not all buffers have Site Scores. 
# Impute the missing values using other confounders and binary treatment using mice
 # Keep only needed columns (plus ID/cluster for merging later)
predictors <- c(
  "Z",
  "Site_Score",
  "population_density",
  "percent_hispanic",
  "percent_black",
  "percent_indigenous",
  "percent_asian",
  "percent_renter_occupied",
  "median_household_income",
  "median_house_value",
  "percent_poverty",
  "percent_high_school_grad",
  "median_year_structure_built",
  "metal", 
  "voc",
  "gas",
  "solid",
  "water"
)
keep_cols <- unique(c(
  "Site_EPA_ID",
  predictors
))
buffers_imp <- st_drop_geometry(buffers[, intersect(names(buffers), keep_cols)])
#buffers_imp[,4:ncol(buffers_imp)] <- scale(buffers_imp[,4:ncol(buffers_imp)]) # standardize predictors

meth <- make.method(buffers_imp)
meth[] <- ""              # default: do not impute
meth[] <- ""  # default to no imputation
# Continuous
meth["median_household_income"] <- "pmm" # 1 missing
meth["median_house_value"] <- "pmm" # 2 missing
meth["median_year_structure_built"] <- "pmm" # just 1 missing
binary_vars <- c("metal", "voc", "gas", "solid", "water") # 188 missing of each
meth[binary_vars[binary_vars %in% names(buffers_imp)]] <- "logreg"

# Predictor matrix
pred <- make.predictorMatrix(buffers_imp)
pred[,] <- 0

# Allow each imputed variable to use all predictors
for (v in names(meth)[meth != ""]) {
  pred[v, predictors[predictors %in% names(buffers_imp)]] <- 1
}

# Missingness is modest. Let's just try single imputation? 
set.seed(123)
imp <- mice(
  data = buffers_imp,
  m = 1,              # single imputation
  maxit = 15,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE
)

buffers_completed <- complete(imp, 1)
# Check order is the same
mean(buffers_completed$Site_EPA_ID == buffers$Site_EPA_ID)
# impute the original missing values of buffers with the imputed values
buffers$metal <- buffers_completed$metal
buffers$voc <- buffers_completed$voc
buffers$gas <- buffers_completed$gas
buffers$solid <- buffers_completed$solid
buffers$water <- buffers_completed$water
buffers$median_household_income <- buffers_completed$median_household_income
buffers$median_house_value <- buffers_completed$median_house_value
buffers$median_year_structure_built <- buffers_completed$median_year_structure_built

# 4: Calculate distance matrix and adjacency matrix
dmat <- distm(cbind(buffers$Longitude, buffers$Latitude), fun = distHaversine)
# Convert distances to kilometers
dmat <- dmat / 1000

# Calculate a 5 nearest neighbors matrix using dmat
n <- nrow(dmat)

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

# Create adjacency matrix: two points are adjacent if one is the other's five nearest neighbors
adjacency_matrix <- nn_matrix 

# Convert metal, voc, gas, solid, water back to numeric 
buffers$metal <- as.numeric(as.character(buffers$metal))
buffers$voc <- as.numeric(as.character(buffers$voc))
buffers$gas <- as.numeric(as.character(buffers$gas))
buffers$solid <- as.numeric(as.character(buffers$solid))
buffers$water <- as.numeric(as.character(buffers$water))

# Save design matrix, distance matrix, and adjacency matrix
X <- buffers[c('Site_Score',
               'metal',
               'voc',
               'gas',
               'solid',
               'water',
               'population_density',
               'percent_hispanic',
               'percent_black',
               'percent_indigenous',
               'percent_asian',
               'percent_renter_occupied',
               'median_household_income',
               'median_house_value',
               'percent_poverty',
               'percent_high_school_grad',
               'median_year_structure_built')]
X <- cbind.data.frame(Intercept = 1, X) # Add intercept
# Drop geometry from X
X <- st_drop_geometry(X)
X <- X[,-ncol(X)]
Z <- buffers$Z
clusters <- buffers$cluster
# Save preprocessed data. 
save(X, Z, dmat, clusters, adjacency_matrix, buffers, data, 
     file = 'data/preprocessed_superfunds.RData')
