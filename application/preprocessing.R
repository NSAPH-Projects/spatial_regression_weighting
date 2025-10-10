library(tidycensus)
library(dplyr)
library(sf)
library(arcpullr)
library(matrixStats)
library(geosphere)
library(Matrix)
library(mice)

setwd('../') # Work in parent directory

################## 1: Exposure Data ###################

# As of June 6, 2024, there were 1,340 Superfund sites in the National Priorities List in the United States.
# Thirty-nine additional sites have been proposed for entry on the list, 
# and 457 sites have been cleaned up and removed from the list

# Read in csv from arcgis data table in url
# url <- 'https://epa.maps.arcgis.com/home/item.html?id=c2b7cdff579c41bbba4898400aa38815'
# https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=33cebcdfdd1b4c3a8b51d416956c41f1
url <- 'https://services.arcgis.com/cJ9YHowT8TU7DUyn/arcgis/rest/services/Superfund_National_Priorities_List_(NPL)_Sites_with_Status_Information/FeatureServer/0'
data <- get_spatial_layer(url)

data <- data %>% 
  dplyr::select(Site_Name, 
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
# Subset those whose Deletion Date is NA or after 2001
# Convert Deletion Date to date
data$Deletion_Date <- as.Date(data$Deletion_Date, format = '%m/%d/%Y')
# Subset to the 1611 superfund sites that were not cleaned up before 2001
data <- data %>% 
  filter(Deletion_Date >= '2001-01-01' | is.na(Deletion_Date))
# Treated units: superfund sites that have been cleaned up between 2001 and 2015
data$Z <- ifelse(data$Deletion_Date <= '2015-12-31' & data$Deletion_Date >= '2001-01-01', 1, 0)
# Turn NAs in Z to 0 (not yet treated)
data$Z[is.na(data$Z)] <- 0
summary(data$Z) # 6% of superfund sites are treated

data_sf <- st_as_sf(data, coords = c('Longitude', 'Latitude'), crs = 4326)

# Perform the spatial join
data_sf$cluster <- substr(data_sf$Site_EPA_ID, 1, 2)
data_sf$cluster <- factor(data_sf$cluster)
data_sf$cluster <- as.numeric(data_sf$cluster)
data_sf <- data_sf[order(data_sf$cluster),] # IMPORTANT TO ORDER! 

sf_points_projected <- st_transform(data_sf, crs = 5070)

# Create buffers with a radius of 2 kilometers (2000 meters)
buffers <- st_buffer(sf_points_projected, dist = 2000)

buffers_original_crs <- st_transform(buffers, crs = st_crs(data_sf))

st_write(buffers_original_crs, "shapefiles/buffers_2k.shp")

#################### 2: Confounder data ###################

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

#################################### Merge with buffers  ####################################

# Read in 2k superfund buffers
buffers <- st_read('shapefiles/buffers_2k.shp')

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
    weighted_population = population_weight * total_population # Estimated pop if population was uniformly distributed
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
    total_population = sum(weighted_population, na.rm = TRUE), # Total estimated population if uniform
    buffer_area = first(buffer_area), # Precomputed buffer area
    population_density = total_population / buffer_area, # Population density
    across(
      c(
        'percent_hispanic', 'percent_black', 'percent_white', 'percent_indigenous',
        'percent_asian', 'percent_poverty', 'percent_high_school_grad'
      ),
      ~ sum(.x, na.rm = TRUE) / sum(intersect_area, na.rm = TRUE) # Area-weighted averages (previously weighted above)
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
    )
  )


# View the resulting buffer-level averages
summary(buffer_averages)

# Now n = 1584. Which sites got dropped?
buffers$S_EPA_I[which(!buffers$S_EPA_I %in% buffer_averages$S_EPA_I)] # PR, Guam, Virgin Islands

# Convert population density (ppl/m^2) to ppl./mi^2
buffer_averages$population_density <- buffer_averages$population_density * 2589988.11
save(buffer_averages, file = 'data/buffer_averages.RData')

#################### 3: Merge exposure and confounder data ####################
# rename the column Site_EPA_ID in buffer_averages to S_EPA_I
buffer_averages <- buffer_averages %>%
  rename(S_EPA_I = Site_EPA_ID)
buffers <- buffers %>%
  left_join(
    st_drop_geometry(buffer_averages),
    by = c("S_EPA_I")
  )
buffers <- buffers %>% dplyr::select(S_EPA_I,
                                     Sit_Scr,
                                      Z,
                                      cluster,
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
# Drop the rows with NA census confounder values (PR, Guam, Virgin Islands)
buffers <- buffers[!is.na(buffers$percent_hispanic),] # 1584
# That line wasn't specific to missingness in Hispanic but census confounders more generally
# Order by cluster
buffers <- buffers[order(buffers$cluster),]
# Not all buffers have Site Scores. 
# Impute the missing values using other confounders and binary treatment using mice
 # Keep only needed columns (plus ID/cluster for merging later)
predictors <- c(
  "Z",
  "population_density",
  "percent_hispanic",
  "percent_black",
  "percent_indigenous",
  "percent_asian",
  "median_household_income",
  "median_house_value",
  "percent_poverty",
  "percent_high_school_grad",
  "median_year_built"
)
keep_cols <- unique(c(
  "S_EPA_I",
  "Sit_Scr",
  predictors
))
buffers_imp <- st_drop_geometry(buffers[, intersect(names(buffers), keep_cols)])
buffers_imp[,4:ncol(buffers_imp)] <- scale(buffers_imp[,4:ncol(buffers_imp)]) # standardize predictors

meth <- make.method(buffers_imp)
meth[] <- ""              # default: do not impute
meth['Sit_Scr'] <- "pmm"     # impute only Sit_Scr

pred <- matrix(0, nrow = ncol(buffers_imp), ncol = ncol(buffers_imp),
               dimnames = list(names(buffers_imp), names(buffers_imp)))
pred['Sit_Scr', predictors[predictors %in% names(buffers_imp)]] <- 1

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
mean(buffers_completed$S_EPA_I == buffers$S_EPA_I)
# impute the original missing values of buffers with the imputed values
buffers$Sit_Scr <- buffers_completed$Sit_Scr

buffers <- na.omit(buffers) # one median house value missing
# n = 1583

# 4: Calculate distance matrix and adjacency matrix
dmat <- distm(cbind(buffers$Longitd, buffers$Latitud), fun = distHaversine)
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

# Save design matrix, distance matrix, and adjacency matrix
X <- buffers[c('Sit_Scr',
               'population_density',
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
# Drop geometry from X
X <- st_drop_geometry(X)
X <- X[,-13]
Z <- buffers$Z
clusters <- buffers$cluster
# Save preprocessed data. 
save(X, Z, dmat, clusters, adjacency_matrix, buffers, file = 'data/preprocessed_superfunds.RData')
