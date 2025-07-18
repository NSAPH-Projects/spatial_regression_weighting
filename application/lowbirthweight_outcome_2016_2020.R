library(pdftools)
library(stringr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(readxl)

################################################ Low birth weight in New Jersey ################################################
# 1: Low Birth weight
# https://www-doh.nj.gov/doh-shad/query/builder/birth/BW/LBW.html
bw <- read_xlsx('../data/nj_2016_2020_lowbwt.xlsx') 
# Select Mother's Municipality, Number of Low Birthweight, and Number of Live Births
bw <- bw[,c(2,3,7,8)]
colnames(bw) <- c('Municipality', 'Year', 'Low_Birthweight', 'Num_Births_Lbw')

# 2: Average Birthweight
# https://www-doh.nj.gov/doh-shad/query/builder/birth/BirthBirthCnty/AvgBirthWt.html
avgbw <- read_xlsx('../data/nj_2016_2020_avg_bw.xlsx')
avgbw <- avgbw[,c(2,3,7,8)]
colnames(avgbw) <- c('Municipality', 'Year', 'Total_Birthweight', 'Num_Births_Bwt')

# 3: Infant Mortality
# https://www-doh.nj.gov/doh-shad/query/builder/infantfetal/Infant/InfMortRate.html
inf <- read_xlsx('../data/nj_2016_2020_infant_mortality.xlsx')
inf <- inf[,c(2,3,7,8)]
colnames(inf) <- c('Municipality', 'Year', 'Infant_Mortality', 'Num_Births_Inf')

# 4: Preterm Birth
# https://www-doh.nj.gov/doh-shad/query/builder/birth/Preterm/All.html
pre <- read_xlsx('../data/nj_2016_2020_preterm.xlsx')
pre <- pre[,c(2,3,7,8)]
colnames(pre) <- c('Municipality', 'Year', 'Preterm', 'Num_Births_Pre')

# 5: Average gestational age
# https://www-doh.nj.gov/doh-shad/query/builder/birth/BirthBirthCnty/AvgGestAge.html
avg_gest <- read_xlsx('../data/nj_2016_2020_avg_gest_age.xlsx')
avg_gest <- avg_gest[,c(2,3,7,8)]
colnames(avg_gest) <- c('Municipality', 'Year', 'Avg_Gestational_Age', 'Num_Births_Avg')

# 1: Birthweight: Remove commas from numbers
bw$Low_Birthweight <- as.numeric(str_remove_all(bw$Low_Birthweight, ','))
bw$Num_Births_Lbw <- as.numeric(str_remove_all(bw$Num_Births_Lbw, ','))
# Convert both to numeric
bw$Low_Birthweight <- as.numeric(bw$Low_Birthweight)
bw$Num_Births_Lbw <- as.numeric(bw$Num_Births_Lbw)
# Aggregate data by municipality
bw <- bw %>%
  group_by(Municipality) %>%
  summarise(Low_Birthweight = sum(Low_Birthweight, na.rm = TRUE),
            Num_Births_Lbw = sum(Num_Births_Lbw, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
bw <- bw[!grepl('unknown', bw$Municipality, ignore.case = TRUE),]
# Just get counties for municipalities that exist in multiple counties
bw$County <- ifelse(grepl('\\(.*\\)', bw$Municipality),
                    str_extract(bw$Municipality, '\\(.*\\)'),
                    NA)
# Remove parentheses from County
bw$County <- str_remove_all(bw$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
bw$Municipality <- str_remove_all(bw$Municipality, '\\(.*\\)')
bw$Municipality <- str_trim(bw$Municipality)

# 2: Average Birthweight: Remove commas from numbers
avgbw$Total_Birthweight <- as.numeric(str_remove_all(avgbw$Total_Birthweight, ','))
avgbw$Num_Births_Bwt <- as.numeric(str_remove_all(avgbw$Num_Births_Bwt, ','))
# Convert both to numeric
avgbw$Total_Birthweight <- as.numeric(avgbw$Total_Birthweight)
avgbw$Num_Births_Bwt <- as.numeric(avgbw$Num_Births_Bwt)
# Aggregate data by municipality
avgbw <- avgbw %>%
  group_by(Municipality) %>%
  summarise(Total_Birthweight = sum(Total_Birthweight, na.rm = TRUE),
            Num_Births_Bwt = sum(Num_Births_Bwt, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
avgbw <- avgbw[!grepl('unknown', avgbw$Municipality, ignore.case = TRUE),]
# Just get counties for municipalities that exist in multiple counties
avgbw$County <- ifelse(grepl('\\(.*\\)', avgbw$Municipality),
                       str_extract(avgbw$Municipality, '\\(.*\\)'),
                       NA)
# Remove parentheses from County
avgbw$County <- str_remove_all(avgbw$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
avgbw$Municipality <- str_remove_all(avgbw$Municipality, '\\(.*\\)')
avgbw$Municipality <- str_trim(avgbw$Municipality)

# 3: Infant Mortality: Remove commas from numbers
inf$Infant_Mortality <- as.numeric(str_remove_all(inf$Infant_Mortality, ','))
inf$Num_Births_Inf <- as.numeric(str_remove_all(inf$Num_Births_Inf, ','))
# Convert both to numeric
inf$Infant_Mortality <- as.numeric(inf$Infant_Mortality)
inf$Num_Births_Inf <- as.numeric(inf$Num_Births_Inf)
# Aggregate data by municipality
inf <- inf %>%
  group_by(Municipality) %>%
  summarise(Infant_Mortality = sum(Infant_Mortality, na.rm = TRUE),
            Num_Births_Inf = sum(Num_Births_Inf, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
inf <- inf[!grepl('unknown', inf$Municipality, ignore.case = TRUE),]
# Just get counties for municipalities that exist in multiple counties
inf$County <- ifelse(grepl('\\(.*\\)', inf$Municipality),
                      str_extract(inf$Municipality, '\\(.*\\)'),
                      NA)
# Remove parentheses from County
inf$County <- str_remove_all(inf$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
inf$Municipality <- str_remove_all(inf$Municipality, '\\(.*\\)')
inf$Municipality <- str_trim(inf$Municipality)

# 4: Preterm: Remove commas from numbers
pre$Preterm <- as.numeric(str_remove_all(pre$Preterm, ','))
pre$Num_Births_Pre <- as.numeric(str_remove_all(pre$Num_Births_Pre, ','))
# Convert both to numeric
pre$Preterm <- as.numeric(pre$Preterm)
pre$Num_Births_Pre <- as.numeric(pre$Num_Births_Pre)
# Aggregate data by municipality
pre <- pre %>%
  group_by(Municipality) %>%
  summarise(Preterm = sum(Preterm, na.rm = TRUE),
            Num_Births_Pre = sum(Num_Births_Pre, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
pre <- pre[!grepl('unknown', pre$Municipality, ignore.case = TRUE),]
# Just get counties for municipalities that exist in multiple counties
pre$County <- ifelse(grepl('\\(.*\\)', pre$Municipality),
                      str_extract(pre$Municipality, '\\(.*\\)'),
                      NA)
# Remove parentheses from County
pre$County <- str_remove_all(pre$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
pre$Municipality <- str_remove_all(pre$Municipality, '\\(.*\\)')
pre$Municipality <- str_trim(pre$Municipality)

# 5: Average gestational age: Remove commas from numbers
avg_gest$Avg_Gestational_Age <- as.numeric(str_remove_all(avg_gest$Avg_Gestational_Age, ','))
avg_gest$Num_Births_Avg <- as.numeric(str_remove_all(avg_gest$Num_Births_Avg, ','))
# Convert both to numeric
avg_gest$Avg_Gestational_Age <- as.numeric(avg_gest$Avg_Gestational_Age)
avg_gest$Num_Births_Avg <- as.numeric(avg_gest$Num_Births_Avg)
# Aggregate data by municipality
avg_gest <- avg_gest %>%
  group_by(Municipality) %>%
  summarise(Avg_Gestational_Age = sum(Avg_Gestational_Age, na.rm = TRUE),
            Num_Births_Avg = sum(Num_Births_Avg, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
avg_gest <- avg_gest[!grepl('unknown', avg_gest$Municipality, ignore.case = TRUE),]
# Just get counties for municipalities that exist in multiple counties
avg_gest$County <- ifelse(grepl('\\(.*\\)', avg_gest$Municipality),
                      str_extract(avg_gest$Municipality, '\\(.*\\)'),
                      NA)
# Remove parentheses from County
avg_gest$County <- str_remove_all(avg_gest$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
avg_gest$Municipality <- str_remove_all(avg_gest$Municipality, '\\(.*\\)')
avg_gest$Municipality <- str_trim(avg_gest$Municipality)

birth_outcomes <- full_join(bw, avgbw, by = c('Municipality', 'County'))
birth_outcomes <- full_join(birth_outcomes, inf, by = c('Municipality', 'County'))
birth_outcomes <- full_join(birth_outcomes, pre, by = c('Municipality', 'County'))
birth_outcomes <- full_join(birth_outcomes, avg_gest, by = c('Municipality', 'County'))

# Manually redo some municipality names so can merge with NJ municipality shapefile
# Convert Lyndhurst Borough to Township
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Lyndhurst Borough',
                                     'Lyndhurst Township',
                                     birth_outcomes$Municipality)
# Convert Caldwell Boro Township to Caldwell Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Caldwell Boro Township',
                                     'Caldwell Borough',
                                     birth_outcomes$Municipality)
# Convert Essex Fells Township to Essex Fells Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Essex Fells Township',
                                     'Essex Fells Borough',
                                     birth_outcomes$Municipality)
# Convert Glen Ridge Boro Township to Glen Ridge Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Glen Ridge Boro Township',
                                     'Glen Ridge Borough',
                                     birth_outcomes$Municipality)
# Convert Orange City to City of Orange Township
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Orange City',
                                     'City of Orange Township',
                                     birth_outcomes$Municipality)
# Convert Princeton Borough to Princeton
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Princeton Borough',
                                     'Princeton',
                                     birth_outcomes$Municipality)
# Convert Lavalette Borough to Lavallette Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Lavalette Borough',
                                     'Lavallette Borough',
                                     birth_outcomes$Municipality)
# Convert Peapack and Gladstone Borough to Peapack-Gladstone Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Peapack and Gladstone Borough',
                                     'Peapack-Gladstone Borough',
                                     birth_outcomes$Municipality)
# Convert Avon-By-The-Sea Borough to Avon-by-the-Sea Borough
birth_outcomes$Municipality <- ifelse(birth_outcomes$Municipality == 'Avon-By-The-Sea Borough',
                                     'Avon-by-the-Sea Borough',
                                     birth_outcomes$Municipality)

# Read in NJ Municipality shapefiles
# https://njogis-newjersey.opendata.arcgis.com/datasets/municipal-boundaries-of-nj-hosted-3424/explore
nj_mun <- st_read('../shapefiles/NJ_Municipal_Boundaries/NJ_Municipal_Boundaries_3424.shp')
mean(birth_outcomes$Municipality %in% nj_mun$MUN_LABEL) # 0.986
mean(nj_mun$MUN_LABEL %in% birth_outcomes$Municipality) # 1

# What are the column names of birth_outcomes minus Municipality and County?
colnames(birth_outcomes)[colnames(birth_outcomes) %in% c('Municipality', 'County') == FALSE]
# Make these new columns NA in nj_mun
nj_mun[, colnames(birth_outcomes)[colnames(birth_outcomes) %in% c('Municipality', 'County') == FALSE]] <- NA

for (i in 1:nrow(birth_outcomes)){
  if (birth_outcomes$Municipality[i] %in% nj_mun$MUN_LABEL){
    # If one match for municipality name, 
    if (is.na(birth_outcomes$County[i])){
      j <- which(nj_mun$MUN_LABEL == birth_outcomes$Municipality[i])
      nj_mun$Low_Birthweight[j] <- birth_outcomes$Low_Birthweight[i]
      nj_mun$Num_Births_Lbw[j] <- birth_outcomes$Num_Births_Lbw[i]
      nj_mun$Total_Birthweight[j] <- birth_outcomes$Total_Birthweight[i]
      nj_mun$Num_Births_Bwt[j] <- birth_outcomes$Num_Births_Bwt[i]
      nj_mun$Infant_Mortality[j] <- birth_outcomes$Infant_Mortality[i]
      nj_mun$Num_Births_Inf[j] <- birth_outcomes$Num_Births_Inf[i]
      nj_mun$Preterm[j] <- birth_outcomes$Preterm[i]
      nj_mun$Num_Births_Pre[j] <- birth_outcomes$Num_Births_Pre[i]
      nj_mun$Avg_Gestational_Age[j] <- birth_outcomes$Avg_Gestational_Age[i]
      nj_mun$Num_Births_Avg[j] <- birth_outcomes$Num_Births_Avg[i]
    }
    # If two matches for municipality name, choose the one with the same county
    if (!is.na(birth_outcomes$County[i])){
      j <- which(nj_mun$MUN_LABEL == birth_outcomes$Municipality[i] & nj_mun$COUNTY == toupper(birth_outcomes$County[i]))
      nj_mun$Low_Birthweight[j] <- birth_outcomes$Low_Birthweight[i]
      nj_mun$Num_Births_Lbw[j] <- birth_outcomes$Num_Births_Lbw[i]
      nj_mun$Total_Birthweight[j] <- birth_outcomes$Total_Birthweight[i]
      nj_mun$Num_Births_Bwt[j] <- birth_outcomes$Num_Births_Bwt[i]
      nj_mun$Infant_Mortality[j] <- birth_outcomes$Infant_Mortality[i]
      nj_mun$Num_Births_Inf[j] <- birth_outcomes$Num_Births_Inf[i]
      nj_mun$Preterm[j] <- birth_outcomes$Preterm[i]
      nj_mun$Num_Births_Pre[j] <- birth_outcomes$Num_Births_Pre[i]
      nj_mun$Avg_Gestational_Age[j] <- birth_outcomes$Avg_Gestational_Age[i]
      nj_mun$Num_Births_Avg[j] <- birth_outcomes$Num_Births_Avg[i]
    }
  }
}

# Plot low birth weight rate on map using geom_sf
glbw <- ggplot(nj_mun) +
  geom_sf(aes(fill = 100*nj_mun$Low_Birthweight/nj_mun$Num_Births_Lbw), color = NA) +
  scale_fill_viridis_c(name = "%") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Plot average birth weight
gavgbw <- ggplot(nj_mun) +
  geom_sf(aes(fill = (1/453.592)*nj_mun$Total_Birthweight/nj_mun$Num_Births_Bwt), color = NA) +
  scale_fill_viridis_c(name = "Pounds") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Average Birth weight") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Plot infant mortality
ginf <- ggplot(nj_mun) +
  geom_sf(aes(fill = 100*nj_mun$Infant_Mortality/nj_mun$Num_Births_Inf), color = NA) +
  scale_fill_viridis_c(name = "%") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Infant Mortality Rate") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Plot preterm birth
gpre <- ggplot(nj_mun) +
  geom_sf(aes(fill = 100*nj_mun$Preterm/nj_mun$Num_Births_Pre), color = NA) +
  scale_fill_viridis_c(name = "%") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Preterm Birth Rate") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Plot gestational age
ggest <- ggplot(nj_mun) +
  geom_sf(aes(fill = nj_mun$Avg_Gestational_Age/nj_mun$Num_Births_Avg), color = NA) +
  scale_fill_viridis_c(name = "weeks") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL) +
  labs(title = "Average Gestational Age") + # remove gridlines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

png('images/birth_outcomes_map.png', width = 1400, height = 1000, res = 150)
grid.arrange(glbw, gavgbw, ginf, gpre, ggest, ncol = 3)
dev.off()

nj_mun <- select(nj_mun, COUNTY, MUN_LABEL, 
                 Low_Birthweight, Num_Births_Lbw,
                 Total_Birthweight, Num_Births_Bwt,
                 Infant_Mortality, Num_Births_Inf,
                 Preterm, Num_Births_Pre,
                 Avg_Gestational_Age, Num_Births_Avg)

colnames(nj_mun)[2] <- 'Municipality'
mun <- nj_mun

st_write(mun, '../shapefiles/Municipality_Birth_Statistics_2016_2020.shp')
