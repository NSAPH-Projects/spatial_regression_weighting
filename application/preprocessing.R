library(tidycensus)
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
library(purrr)
library(stringr)
library(janitor)
library(lubridate)
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
#################### 2: Confounder data ###################

# # Read in fips_codes
# data(fips_codes)
# statefps <- unique(fips_codes$state_code)
# statenames <- unique(fips_codes$state_name)
# fips_codes <- cbind.data.frame(state_code = statefps, state_name = statenames)
# # remove the last 6 rows (not states)
# fips_codes <- fips_codes[-c((nrow(fips_codes)-5):nrow(fips_codes)),]
# 
# setwd('shapefiles/')
# all_shapefiles <- list()
# # Unzip and read in each state's tract shapefile in 2000 that has the form: https://www2.census.gov/geo/pvs/tiger2010st/08_Colorado/08/tl_2010_08_tract00.zip
# baseURL <- "https://www2.census.gov/geo/pvs/tiger2010st/"
# for (state in fips_codes$state_code) {
#   state_name <- fips_codes$state_name[fips_codes$state_code == state]
#   state_name <- gsub(" ", "_", state_name)
#   state_name <- gsub(",", "", state_name)
#   state_name <- gsub(" ", "", state_name)
#   state_name <- paste0(state, "_", state_name)
#   state_name <- paste0(state_name, "/", state, "/")
#   url <- paste0(baseURL, state_name, "tl_2010_", state, "_tract00.zip")
#   print(url)
#   download.file(url, destfile = paste0(state, "_tract00.zip"))
#   unzip(paste0(state, "_tract00.zip"))
#   # Read in shapefile from the unzipped folder
#   state_shapefile <- st_read(paste0('tl_2010_', state, "_tract00.shp"))
#   all_shapefiles[[state]] <- state_shapefile
# }
# # Combine all shapefiles into one
# all_shapefiles <- do.call(rbind, all_shapefiles)
# setwd('../')
# 
# # Define the list of variables for the 2000 Census
# vars_2000_sf1 <- c(
#   total_population = "P001001",           # Total population
#   race_denom = "P009001",                 # Total population for race data
#   white_population = "P009002",           # White population"
#   black_population = "P009003",           # Black population
#   indigenous_population = "P009004",      # American Indian and Alaska Native population
#   asian_population = "P009005",           # Asian population
#   hispanic_denom = "P004001",             # Total population for Hispanic data
#   hispanic_population = "P004002"        # Hispanic population
# )
# 
# vars_2000_sf3 <- c(
#   median_household_income = "P053001",    # Median household income in 1999
#   median_house_value = "H085001",         # Median house value
#   poverty_population = "P087002",         # Population with income below poverty level in 1999
#   poverty_denom = "P087001",              # Total population for poverty rate calculation
#   high_school_grad_male = "P037011",      # Male Population with HS degree
#   some_college_1_male = "P037012",        # Male Population with <1 year of college
#   some_college_male = "P037013",          # Male Population with some college
#   college_association_male = "P037014",   # Male Population with Associate's degree
#   college_bachelor_male = "P037015",      # Male Population with Bachelor's degree
#   college_masters_male = "P037016",       # Male Population with Master's degree
#   college_prof_male = "P037017",          # Male Population with Professional degree
#   college_doctorate_male = "P037018",     # Male Population with Doctorate degree
#   high_school_grad_female = "P037028",    # Female Population with HS degree
#   some_college_1_female = "P037029",      # Female Population with <1 year of college
#   some_college_female = "P037030",        # Female Population with some college
#   college_association_female = "P037031",   # Female Population with Associate's degree
#   college_bachelor_female = "P037032",      # Female Population with Bachelor's degree
#   college_masters_female = "P037033",       # Female Population with Master's degree
#   college_prof_female = "P037034",          # Female Population with Professional degree
#   college_doctorate_female = "P037035",     # Female Population with Doctorate degree
#   high_school_grad_male_denom = "P037002",     # Total male population aged 25 or older with education info
#   high_school_grad_female_denom = "P037019",   # Total female population aged 25 or older with education info
#   median_year_structure_built = "H035001"           # Median year housing units were built
# )
# 
# # Get tract-level data from the 2000 Census
# # Initialize empty data frame of 10 columns
# tract_data_2000_sf1 <- data.frame(matrix(ncol = 10, nrow = 0))
# # Loop through states
# for (state in statefps) {
#   # Get data for each state
#   state_data <- get_decennial(
#     geography = "tract",
#     year = 2000,
#     variables = vars_2000_sf1,
#     state = state,
#     sumfile = 'sf1',
#     output = "wide"
#   )
#   # Append state data to the empty data frame
#   tract_data_2000_sf1 <- rbind(tract_data_2000_sf1, state_data)
# }
# 
# tract_data_2000_sf3 <- data.frame(matrix(ncol = 10, nrow = 0))
# # Loop through states
# for (state in statefps) {
#   # Get data for each state
#   state_data <- get_decennial(
#     geography = "tract",
#     year = 2000,
#     variables = vars_2000_sf3,
#     state = state,
#     sumfile = 'sf3',
#     output = "wide"
#   )
#   # Append state data to the empty data frame
#   tract_data_2000_sf3 <- rbind(tract_data_2000_sf3, state_data)
# }
# 
# tract_data_2000 <- left_join(tract_data_2000_sf1, 
#                              tract_data_2000_sf3, 
#                              by = c("GEOID", "NAME"))
# 
# # Calculate derived variables (percentages)
# tract_data_2000 <- tract_data_2000 %>%
#   mutate(
#     percent_hispanic = hispanic_population/hispanic_denom,
#     percent_black = black_population/race_denom,
#     percent_white = white_population/race_denom,
#     percent_indigenous = indigenous_population/race_denom,
#     percent_asian = asian_population/race_denom,
#     percent_poverty = poverty_population / poverty_denom,
#     percent_high_school_grad = (
#       high_school_grad_male + high_school_grad_female + some_college_1_male + some_college_1_female + some_college_male + some_college_female 
#       + college_association_male + college_association_female + college_bachelor_male + college_bachelor_female 
#       + college_masters_male + college_masters_female + college_prof_male + college_prof_female 
#       + college_doctorate_male + college_doctorate_female
#     ) / (high_school_grad_male_denom + high_school_grad_female_denom)
#   ) %>%
#   select(
#     GEOID, NAME, total_population, 
#     percent_hispanic, percent_black, percent_white, percent_indigenous,
#     percent_asian,  
#     median_household_income, median_house_value, 
#     percent_poverty, percent_high_school_grad, 
#     median_year_structure_built
#   )
# 
# # View the resulting data
# head(tract_data_2000)
# summary(tract_data_2000)
# # Merge tract_data_2000 with all_shapefiles
# tract_data_2000 <- left_join(all_shapefiles, 
#                              tract_data_2000, 
#                              by = c("CTIDFP00" = "GEOID"))
# # change nonsensical values to NA
# tract_data_2000 <- tract_data_2000 %>% mutate(
#   median_house_value = ifelse(median_house_value == 0, NA, median_house_value),
#   median_household_income = ifelse(median_household_income == 0, NA, median_household_income),
#   median_year_structure_built = ifelse(median_year_structure_built == 0, NA, median_year_structure_built)
# )
# save(tract_data_2000, file = 'data/tract_data_2000.RData') # 65443 tracts. 

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
