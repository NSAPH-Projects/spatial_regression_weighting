# Load necessary libraries
library(pdftools)
library(stringr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(readxl)

# Define the file path
# Source: https://www.pa.gov/content/dam/copapwp-pagov/en/health/documents/topics/healthstatistics/vitalstatistics/birthstatistics/documents/Birth_BirthWt_MCD_2016_2020.pdf
file_path <- "../data/PA_Birth_BirthWt_MCD_2016_2020.pdf"

# Extract text from the PDF
pdf_text <- pdf_text(file_path)

# Initialize an empty dataframe to store extracted birth data
birth_data <- data.frame(
  County = character(),
  Municipality = character(),
  Total_Births = integer(),
  Less_1500_Grams = integer(),
  Grams_1500_2499 = integer(),
  Grams_2500_plus = integer(),
  Unknown = integer(),
  stringsAsFactors = FALSE
)

# Define the pattern to detect county names (assuming they appear alone and in uppercase)
county_pattern <- "^[A-Z][A-Za-z\\s]+$"

# Process each page of the PDF
for (page_num in seq_along(pdf_text)) {
  page <- pdf_text[[page_num]]
  lines <- str_split(page, "\n")[[1]]
  county <- NA  # Reset county for each page
  
  for (line_num in seq_along(lines)) {
    line <- str_trim(lines[[line_num]])
    
    # Debugging: Print detected county names
    if (grepl(county_pattern, line) && !grepl("Total", line)) {
      print(paste("Detected County:", line, "on page", page_num))
      county <- line  # Update the current county
      next  # Move to the next line
    }
    
    # Municipality data lines: Assume format "Municipality Name 1,234 12 34 5,678 7"
    if (!is.na(county) && grepl("^\\D+\\s+\\d+", line)) {
      # Remove commas from numbers
      line <- gsub(",", "", line)
      
      parts <- str_match(line, "^(.*?)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)$")
      
      if (!is.na(parts[1, 1])) {
        # Successfully extracted data
        birth_data <- rbind(birth_data, data.frame(
          County = county,
          Municipality = str_trim(parts[1, 2]),
          Total_Births = as.integer(parts[1, 3]),
          Less_1500_Grams = as.integer(parts[1, 4]),
          Grams_1500_2499 = as.integer(parts[1, 5]),
          Grams_2500_plus = as.integer(parts[1, 6]),
          Unknown = as.integer(parts[1, 7]),
          stringsAsFactors = FALSE
        ))
      } else {
        # Debugging: Print lines where extraction failed
        print(paste("Failed to extract municipality data on page", page_num, 
                    "line", line_num, ":", line))
      }
    }
  }
}

# Display the first few rows of the extracted data
head(birth_data)

# Add Wheatland Borough in Mercer to Hermitage City in Mercer
birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] = birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] +
  birth_data[which(birth_data$Municipality == "Wheatland Borough"),3:7]
# remove Wheatland Borough
birth_data = subset(birth_data, Municipality != "Wheatland Borough")
summary(birth_data)

# Remove any municipalities that are "Total " + word (e.g. "Total Potter")
birth_data <- birth_data[!grepl("^Total\\s", birth_data$Municipality),]
write.csv(birth_data, "../data/birth_data.csv", row.names = FALSE)

# Drop the rows where Total_Births = 1
birth_data <- read.csv("../data/birth_data.csv")
#birth_data <- birth_data[birth_data$Total_Births > 1,]

# Calculate the low birth weight rate
birth_data$Low_Birthweight <- birth_data$Less_1500_Grams + birth_data$Grams_1500_2499
birth_data$Num_Births <- birth_data$Total_Births - birth_data$Unknown
birth_data <- select(birth_data, County, Municipality, Low_Birthweight, Num_Births)

# Remove any periods or apostrophes
birth_data$Municipality <- gsub("[[:punct:]]", "", birth_data$Municipality)
# Remoe the space in Dubois
birth_data$Municipality <- gsub("Du Bois City", "DuBois City", birth_data$Municipality)
# Add a space in Bradfordwoods
birth_data$Municipality <- gsub("Bradfordwoods Borough", "Bradford Woods Borough", birth_data$Municipality)
# Remove space in Black Lick
birth_data$Municipality <- gsub("Black Lick Township", "BlackLick Township", birth_data$Municipality)
# Add space to Nantyglo
birth_data$Municipality <- gsub("Nantyglo Borough", "Nanty Glo Borough", birth_data$Municipality)
# Add space to Oilcreek
birth_data$Municipality <- gsub("Oilcreek Township", "Oil Creek Township", birth_data$Municipality)
# Add hyphen to Wilkesbarre
birth_data$Municipality <- gsub("Wilkesbarre City", "Wilkes-Barre City", birth_data$Municipality)
# Add hyphen to Wilkesbarre
birth_data$Municipality <- gsub("Wilkesbarre Township", "Wilkes-Barre Township", birth_data$Municipality)
# Change St Clair Borough to to Saint Clair Borough
birth_data$Municipality <- gsub("St Clair Borough", "Saint Clair Borough", birth_data$Municipality)

# Merge with PA shapefile.
# Source: https://www.pasda.psu.edu/uci/DataSummary.aspx?dataset=41
pamun <- st_read('../shapefiles/PaMunicipalities2025_01.shp')
# to each municipal1 name, add township for CLASS_OF_M = 1TWP 2TWP, borough for BORO, City for CITY, town for TOWN 
pamun$MUNICIPAL1 <- ifelse((pamun$CLASS_OF_M == "1TWP" | pamun$CLASS_OF_M == "2TWP"), 
                           paste0(pamun$MUNICIPAL1, " TOWNSHIP"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "BORO", 
                           paste0(pamun$MUNICIPAL1, " BOROUGH"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "CITY",
                           paste0(pamun$MUNICIPAL1, " CITY"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "TOWN",
                           paste0(pamun$MUNICIPAL1, " TOWN"), pamun$MUNICIPAL1)

mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1)
birth_data[which(!toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1),]

# Some of the non-matches are due to "Mount" in the birth_data Muncipality name but "MT" in the shapefile.
# Fix this by replacing "Mount" with "MT" in the birth_data Municipality name, just for these non-matches.
birth_data$Municipality <- ifelse(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1, 
                                   birth_data$Municipality, 
                                   gsub("Mount", "MT", birth_data$Municipality))
birth_data$Municipality[birth_data$Municipality == "Murrysville Municipality"] = "MURRYSVILLE BOROUGH"
birth_data$Municipality[birth_data$Municipality == "Monroeville Municipality" & birth_data$County == 'Allegheny'] = "MONROEVILLE BOROUGH"
birth_data$Municipality[birth_data$Municipality == "Bethel Park Municipality" & birth_data$County == 'Allegheny'] = "BETHEL PARK BOROUGH"

mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1) # 1. YAY! All match.
birth_data$Municipality <- toupper(birth_data$Municipality)

mean(toupper(birth_data$County) %in% pamun$COUNTY_NAM) # YAY! All match.
birth_data$County <- toupper(birth_data$County)

# Which pamun$MUNICIPAL1 are not in birth_data$Municipality
st_drop_geometry(pamun[which(!(pamun$MUNICIPAL1 %in% birth_data$Municipality)),c('MUNICIPAL1', 'COUNTY_NAM')]) # THERE ARE 14 missing!
# I can't find bw for these.

pamun <- left_join(pamun, birth_data, by = c("COUNTY_NAM" = "County", "MUNICIPAL1" = "Municipality"))

# # Plot low birth weight rate on map using geom_sf
ggplot(pamun) +
  geom_sf(aes(fill = pamun$Low_Birthweight/pamun$Num_Births), color = "black") +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in Pennsylvania")


################################################ Low birth weight in New Jersey ################################################
# Source: https://www-doh.nj.gov/doh-shad/query/builder/birth/BW/LBW.html
bw <- read_xlsx('../data/nj_2016_2020_lowbwt.xlsx') 
# Select Mother's Municipality, Number of Low Birthweight, and Number of Live Births
bw <- bw[,c(2,3,7,8)]
colnames(bw) <- c('Municipality', 'Year', 'Low_Birthweight', 'Num_Births')
# Low Birth weight rate
# Remove commas from numbers
bw$Low_Birthweight <- as.numeric(str_remove_all(bw$Low_Birthweight, ','))
bw$Num_Births <- as.numeric(str_remove_all(bw$Num_Births, ','))
# Convert both to numeric
bw$Low_Birthweight <- as.numeric(bw$Low_Birthweight)
bw$Num_Births <- as.numeric(bw$Num_Births)
# Aggregate data by municipality
bw <- bw %>%
  group_by(Municipality) %>%
  summarise(Low_Birthweight = sum(Low_Birthweight, na.rm = TRUE),
            Num_Births = sum(Num_Births, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
bw <- bw[!grepl('unknown', bw$Municipality, ignore.case = TRUE),]
bw$County <- ifelse(grepl('\\(.*\\)', bw$Municipality),
                    str_extract(bw$Municipality, '\\(.*\\)'),
                    NA)
# Remove parentheses from County
bw$County <- str_remove_all(bw$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
bw$Municipality <- str_remove_all(bw$Municipality, '\\(.*\\)')
bw$Municipality <- str_trim(bw$Municipality)

# Convert Lyndhurst Borough to Township
bw$Municipality <- ifelse(bw$Municipality == 'Lyndhurst Borough',
                                     'Lyndhurst Township',
                                     bw$Municipality)
# Convert Caldwell Boro Township to Caldwell Borough
bw$Municipality <- ifelse(bw$Municipality == 'Caldwell Boro Township',
                                     'Caldwell Borough',
                                     bw$Municipality)
# Convert Essex Fells Township to Essex Fells Borough
bw$Municipality <- ifelse(bw$Municipality == 'Essex Fells Township',
                                     'Essex Fells Borough',
                                     bw$Municipality)
# Convert Glen Ridge Boro Township to Glen Ridge Borough
bw$Municipality <- ifelse(bw$Municipality == 'Glen Ridge Boro Township',
                                     'Glen Ridge Borough',
                                     bw$Municipality)
# Convert Orange City to City of Orange Township
bw$Municipality <- ifelse(bw$Municipality == 'Orange City',
                                     'City of Orange Township',
                                     bw$Municipality)
# Convert Princeton Borough to Princeton
bw$Municipality <- ifelse(bw$Municipality == 'Princeton Borough',
                                     'Princeton',
                                     bw$Municipality)
# Convert Lavalette Borough to Lavallette Borough
bw$Municipality <- ifelse(bw$Municipality == 'Lavalette Borough',
                                     'Lavallette Borough',
                                     bw$Municipality)
# Convert Peapack and Gladstone Borough to Peapack-Gladstone Borough
bw$Municipality <- ifelse(bw$Municipality == 'Peapack and Gladstone Borough',
                                     'Peapack-Gladstone Borough',
                                     bw$Municipality)
# Convert Avon-By-The-Sea Borough to Avon-by-the-Sea Borough
bw$Municipality <- ifelse(bw$Municipality == 'Avon-By-The-Sea Borough',
                                     'Avon-by-the-Sea Borough',
                                     bw$Municipality)

# Read in NJ Municipality shapefiles
# https://njogis-newjersey.opendata.arcgis.com/datasets/municipal-boundaries-of-nj-hosted-3424/explore
nj_mun <- st_read('../shapefiles/NJ_Municipal_Boundaries/NJ_Municipal_Boundaries_3424.shp')
mean(bw$Municipality %in% nj_mun$MUN_LABEL) # 0.986
mean(nj_mun$MUN_LABEL %in% bw$Municipality) # 1

#bw$Low_Birthweight_Rate <- bw$Low_Birthweight / bw$Num_Births
#hist(bw$Low_Birthweight_Rate)

nj_mun$Low_Birthweight <- NA
nj_mun$Num_Births <- NA
for (i in 1:nrow(bw)){
  if (bw$Municipality[i] %in% nj_mun$MUN_LABEL){
    # If one match for municipality name, 
    if (is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$Municipality[i])
      nj_mun$Low_Birthweight[j] <- bw$Low_Birthweight[i]
      nj_mun$Num_Births[j] <- bw$Num_Births[i]
    }
    # If two matches for municipality name, choose the one with the same county
    if (!is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$Municipality[i] & nj_mun$COUNTY == toupper(bw$County[i]))
      nj_mun$Low_Birthweight[j] <- bw$Low_Birthweight[i]
      nj_mun$Num_Births[j] <- bw$Num_Births[i]
    }
  }
}

# Plot low birth weight rate on map using geom_sf
ggplot(nj_mun) +
  geom_sf(aes(fill = nj_mun$Low_Birthweight/nj_mun$Num_Births), color = NA) +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in New Jersey")

# Combine the PA shapefile with the NJ shapefile
pamun <- st_transform(pamun, crs = st_crs(nj_mun))
pamun <- select(pamun, COUNTY, MUNICIPAL1, Low_Birthweight, Num_Births)
nj_mun <- select(nj_mun, COUNTY, MUN_LABEL, Low_Birthweight, Num_Births)
# Rename MUNICIPAL1 and MUN_LABEL each to Municipality
colnames(pamun)[2] <- 'Municipality'
colnames(nj_mun)[2] <- 'Municipality'
# Combine the two dataframes
mun <- rbind(pamun, nj_mun)

# Plot low birth weight rate on map using geom_sf
ggplot(mun) +
  geom_sf(aes(fill = mun$Low_Birthweight/mun$Num_Births), color = NA) +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in New Jersey and Pennsylvania")
st_write(mun, '../shapefiles/Municipality_Low_Birthweight_Rate_2016_2020.shp')

# ################################################ VERY low birth weight in Pennsylvania #############################################
# 
# # Define the file path
# # Source: https://www.pa.gov/content/dam/copapwp-pagov/en/health/documents/topics/healthstatistics/vitalstatistics/birthstatistics/documents/Birth_BirthWt_MCD_2016_2020.pdf
# file_path <- "../data/PA_Birth_BirthWt_MCD_2016_2020.pdf"
# 
# # Extract text from the PDF
# pdf_text <- pdf_text(file_path)
# 
# # Initialize an empty dataframe to store extracted birth data
# birth_data <- data.frame(
#   County = character(),
#   Municipality = character(),
#   Total_Births = integer(),
#   Less_1500_Grams = integer(),
#   Grams_1500_2499 = integer(),
#   Grams_2500_plus = integer(),
#   Unknown = integer(),
#   stringsAsFactors = FALSE
# )
# 
# # Define the pattern to detect county names (assuming they appear alone and in uppercase)
# county_pattern <- "^[A-Z][A-Za-z\\s]+$"
# 
# # Process each page of the PDF
# for (page_num in seq_along(pdf_text)) {
#   page <- pdf_text[[page_num]]
#   lines <- str_split(page, "\n")[[1]]
#   county <- NA  # Reset county for each page
#   
#   for (line_num in seq_along(lines)) {
#     line <- str_trim(lines[[line_num]])
#     
#     # Debugging: Print detected county names
#     if (grepl(county_pattern, line) && !grepl("Total", line)) {
#       print(paste("Detected County:", line, "on page", page_num))
#       county <- line  # Update the current county
#       next  # Move to the next line
#     }
#     
#     # Municipality data lines: Assume format "Municipality Name 1,234 12 34 5,678 7"
#     if (!is.na(county) && grepl("^\\D+\\s+\\d+", line)) {
#       # Remove commas from numbers
#       line <- gsub(",", "", line)
#       
#       parts <- str_match(line, "^(.*?)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)$")
#       
#       if (!is.na(parts[1, 1])) {
#         # Successfully extracted data
#         birth_data <- rbind(birth_data, data.frame(
#           County = county,
#           Municipality = str_trim(parts[1, 2]),
#           Total_Births = as.integer(parts[1, 3]),
#           Less_1500_Grams = as.integer(parts[1, 4]),
#           Grams_1500_2499 = as.integer(parts[1, 5]),
#           Grams_2500_plus = as.integer(parts[1, 6]),
#           Unknown = as.integer(parts[1, 7]),
#           stringsAsFactors = FALSE
#         ))
#       } else {
#         # Debugging: Print lines where extraction failed
#         print(paste("Failed to extract municipality data on page", page_num, 
#                     "line", line_num, ":", line))
#       }
#     }
#   }
# }
# 
# # Display the first few rows of the extracted data
# head(birth_data)
# 
# # Add Wheatland Borough in Mercer to Hermitage City in Mercer
# birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] = birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] +
#   birth_data[which(birth_data$Municipality == "Wheatland Borough"),3:7]
# # remove Wheatland Borough
# birth_data = subset(birth_data, Municipality != "Wheatland Borough")
# summary(birth_data)
# 
# # # Remove any municipalities that are "Total " + word (e.g. "Total Potter")
# # birth_data <- birth_data[!grepl("^Total\\s", birth_data$Municipality),]
# # #write.csv(birth_data, "../data/birth_data.csv", row.names = FALSE)
# # 
# # # Drop the rows where Total_Births = 1
# # birth_data <- read.csv("../data/birth_data.csv")
# # #birth_data <- birth_data[birth_data$Total_Births > 1,]
# 
# # Calculate the low birth weight rate
# birth_data$VLow_Birthweight <- birth_data$Less_1500_Grams 
# birth_data$Num_Births <- birth_data$Total_Births - birth_data$Unknown
# birth_data <- select(birth_data, County, Municipality, VLow_Birthweight, Num_Births)
# 
# # Remove any periods or apostrophes
# birth_data$Municipality <- gsub("[[:punct:]]", "", birth_data$Municipality)
# # Remoe the space in Dubois
# birth_data$Municipality <- gsub("Du Bois City", "DuBois City", birth_data$Municipality)
# # Add a space in Bradfordwoods
# birth_data$Municipality <- gsub("Bradfordwoods Borough", "Bradford Woods Borough", birth_data$Municipality)
# # Remove space in Black Lick
# birth_data$Municipality <- gsub("Black Lick Township", "BlackLick Township", birth_data$Municipality)
# # Add space to Nantyglo
# birth_data$Municipality <- gsub("Nantyglo Borough", "Nanty Glo Borough", birth_data$Municipality)
# # Add space to Oilcreek
# birth_data$Municipality <- gsub("Oilcreek Township", "Oil Creek Township", birth_data$Municipality)
# # Add hyphen to Wilkesbarre
# birth_data$Municipality <- gsub("Wilkesbarre City", "Wilkes-Barre City", birth_data$Municipality)
# # Add hyphen to Wilkesbarre
# birth_data$Municipality <- gsub("Wilkesbarre Township", "Wilkes-Barre Township", birth_data$Municipality)
# # Change St Clair Borough to to Saint Clair Borough
# birth_data$Municipality <- gsub("St Clair Borough", "Saint Clair Borough", birth_data$Municipality)
# 
# # Merge with PA shapefile.
# # Source: https://www.pasda.psu.edu/uci/DataSummary.aspx?dataset=41
# pamun <- st_read('../shapefiles/PaMunicipalities2025_01.shp')
# # to each municipal1 name, add township for CLASS_OF_M = 1TWP 2TWP, borough for BORO, City for CITY, town for TOWN 
# pamun$MUNICIPAL1 <- ifelse((pamun$CLASS_OF_M == "1TWP" | pamun$CLASS_OF_M == "2TWP"), 
#                            paste0(pamun$MUNICIPAL1, " TOWNSHIP"), pamun$MUNICIPAL1)
# pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "BORO", 
#                            paste0(pamun$MUNICIPAL1, " BOROUGH"), pamun$MUNICIPAL1)
# pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "CITY",
#                            paste0(pamun$MUNICIPAL1, " CITY"), pamun$MUNICIPAL1)
# pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "TOWN",
#                            paste0(pamun$MUNICIPAL1, " TOWN"), pamun$MUNICIPAL1)
# 
# mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1)
# birth_data[which(!toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1),]
# 
# # Some of the non-matches are due to "Mount" in the birth_data Muncipality name but "MT" in the shapefile.
# # Fix this by replacing "Mount" with "MT" in the birth_data Municipality name, just for these non-matches.
# birth_data$Municipality <- ifelse(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1, 
#                                   birth_data$Municipality, 
#                                   gsub("Mount", "MT", birth_data$Municipality))
# birth_data$Municipality[birth_data$Municipality == "Murrysville Municipality"] = "MURRYSVILLE BOROUGH"
# birth_data$Municipality[birth_data$Municipality == "Monroeville Municipality" & birth_data$County == 'Allegheny'] = "MONROEVILLE BOROUGH"
# birth_data$Municipality[birth_data$Municipality == "Bethel Park Municipality" & birth_data$County == 'Allegheny'] = "BETHEL PARK BOROUGH"
# 
# mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1) # 1. YAY! All match.
# birth_data$Municipality <- toupper(birth_data$Municipality)
# 
# mean(toupper(birth_data$County) %in% pamun$COUNTY_NAM) # YAY! All match.
# birth_data$County <- toupper(birth_data$County)
# 
# # Which pamun$MUNICIPAL1 are not in birth_data$Municipality
# st_drop_geometry(pamun[which(!(pamun$MUNICIPAL1 %in% birth_data$Municipality)),c('MUNICIPAL1', 'COUNTY_NAM')]) # THERE ARE 14 missing!
# # I can't find bw for these.
# 
# pamun <- left_join(pamun, birth_data, by = c("COUNTY_NAM" = "County", "MUNICIPAL1" = "Municipality"))
# 
# # # Plot low birth weight rate on map using geom_sf
# ggplot(pamun) +
#   geom_sf(aes(fill = pamun$VLow_Birthweight/pamun$Num_Births), color = "black") +
#   scale_fill_viridis_c(name = "VLow Birthweight Rate") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(title = "VLow Birthweight Rate by Municipality in Pennsylvania")

# ################################################ VERY Low birth weight in New Jersey ################################################
# 
# # Source: https://www-doh.nj.gov/doh-shad/query/builder/birth/BW/VLBW.html
# bw <- read_xlsx('../data/nj_2016_2020_verylowbwt.xlsx') 
# # Select Mother's Municipality, Number of Low Birthweight, and Number of Live Births
# bw <- bw[,c(2,3,7,8)]
# colnames(bw) <- c('Municipality', 'Year', 'VLow_Birthweight', 'Num_Births')
# # Low Birth weight rate
# # Remove commas from numbers
# bw$VLow_Birthweight <- as.numeric(str_remove_all(bw$VLow_Birthweight, ','))
# bw$Num_Births <- as.numeric(str_remove_all(bw$Num_Births, ','))
# # Convert both to numeric
# bw$VLow_Birthweight <- as.numeric(bw$VLow_Birthweight)
# bw$Num_Births <- as.numeric(bw$Num_Births)
# # Aggregate data by municipality
# bw <- bw %>%
#   group_by(Municipality) %>%
#   summarise(VLow_Birthweight = sum(VLow_Birthweight, na.rm = TRUE),
#             Num_Births = sum(Num_Births, na.rm = TRUE)) %>%
#   ungroup()
# # If municipality contains the word "unknown", drop the row
# bw <- bw[!grepl('unknown', bw$Municipality, ignore.case = TRUE),]
# bw$County <- ifelse(grepl('\\(.*\\)', bw$Municipality),
#                     str_extract(bw$Municipality, '\\(.*\\)'),
#                     NA)
# # Remove parentheses from County
# bw$County <- str_remove_all(bw$County, '\\(|\\)')
# # Remove the county and the extra space from the municipality name
# bw$Municipality <- str_remove_all(bw$Municipality, '\\(.*\\)')
# bw$Municipality <- str_trim(bw$Municipality)
# 
# # Convert Lyndhurst Borough to Township
# bw$Municipality <- ifelse(bw$Municipality == 'Lyndhurst Borough',
#                           'Lyndhurst Township',
#                           bw$Municipality)
# # Convert Caldwell Boro Township to Caldwell Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Caldwell Boro Township',
#                           'Caldwell Borough',
#                           bw$Municipality)
# # Convert Essex Fells Township to Essex Fells Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Essex Fells Township',
#                           'Essex Fells Borough',
#                           bw$Municipality)
# # Convert Glen Ridge Boro Township to Glen Ridge Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Glen Ridge Boro Township',
#                           'Glen Ridge Borough',
#                           bw$Municipality)
# # Convert Orange City to City of Orange Township
# bw$Municipality <- ifelse(bw$Municipality == 'Orange City',
#                           'City of Orange Township',
#                           bw$Municipality)
# # Convert Princeton Borough to Princeton
# bw$Municipality <- ifelse(bw$Municipality == 'Princeton Borough',
#                           'Princeton',
#                           bw$Municipality)
# # Convert Lavalette Borough to Lavallette Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Lavalette Borough',
#                           'Lavallette Borough',
#                           bw$Municipality)
# # Convert Peapack and Gladstone Borough to Peapack-Gladstone Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Peapack and Gladstone Borough',
#                           'Peapack-Gladstone Borough',
#                           bw$Municipality)
# # Convert Avon-By-The-Sea Borough to Avon-by-the-Sea Borough
# bw$Municipality <- ifelse(bw$Municipality == 'Avon-By-The-Sea Borough',
#                           'Avon-by-the-Sea Borough',
#                           bw$Municipality)
# 
# # Read in NJ Municipality shapefiles
# # https://njogis-newjersey.opendata.arcgis.com/datasets/municipal-boundaries-of-nj-hosted-3424/explore
# nj_mun <- st_read('../shapefiles/NJ_Municipal_Boundaries/NJ_Municipal_Boundaries_3424.shp')
# mean(bw$Municipality %in% nj_mun$MUN_LABEL) # 0.986
# mean(nj_mun$MUN_LABEL %in% bw$Municipality) # 1
# 
# #bw$Low_Birthweight_Rate <- bw$Low_Birthweight / bw$Num_Births
# #hist(bw$Low_Birthweight_Rate)
# 
# nj_mun$VLow_Birthweight <- NA
# nj_mun$Num_Births <- NA
# for (i in 1:nrow(bw)){
#   if (bw$Municipality[i] %in% nj_mun$MUN_LABEL){
#     # If one match for municipality name, 
#     if (is.na(bw$County[i])){
#       j <- which(nj_mun$MUN_LABEL == bw$Municipality[i])
#       nj_mun$VLow_Birthweight[j] <- bw$VLow_Birthweight[i]
#       nj_mun$Num_Births[j] <- bw$Num_Births[i]
#     }
#     # If two matches for municipality name, choose the one with the same county
#     if (!is.na(bw$County[i])){
#       j <- which(nj_mun$MUN_LABEL == bw$Municipality[i] & nj_mun$COUNTY == toupper(bw$County[i]))
#       nj_mun$VLow_Birthweight[j] <- bw$VLow_Birthweight[i]
#       nj_mun$Num_Births[j] <- bw$Num_Births[i]
#     }
#   }
# }
# 
# # Plot low birth weight rate on map using geom_sf
# ggplot(nj_mun) +
#   geom_sf(aes(fill = nj_mun$VLow_Birthweight/nj_mun$Num_Births), color = NA) +
#   scale_fill_viridis_c(name = "VLow Birthweight Rate") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(title = "VLow Birthweight Rate by Municipality in New Jersey")
# 
# 
# # Combine the PA shapefile with the NJ shapefile
# pamun <- st_transform(pamun, crs = st_crs(nj_mun))
# pamun <- select(pamun, COUNTY, MUNICIPAL1, VLow_Birthweight, Num_Births)
# nj_mun <- select(nj_mun, COUNTY, MUN_LABEL, VLow_Birthweight, Num_Births)
# # Rename MUNICIPAL1 and MUN_LABEL each to Municipality
# colnames(pamun)[2] <- 'Municipality'
# colnames(nj_mun)[2] <- 'Municipality'
# # Combine the two dataframes
# mun_vlow <- rbind(pamun, nj_mun)
# 
# # Plot low birth weight rate on map using geom_sf
# ggplot(mun_vlow) +
#   geom_sf(aes(fill = mun_vlow$VLow_Birthweight/mun_vlow$Num_Births), color = NA) +
#   scale_fill_viridis_c(name = "VLow Birthweight Rate") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(title = "VLow Birthweight Rate by Municipality in New Jersey and Pennsylvania")
# st_write(mun_vlow, '../shapefiles/Municipality_VLow_Birthweight_Rate_2016_2020.shp')

######################################## Baseline birth rates (2001-2005) #############################################
file_path <- "../data/PA_Birth_BirthWt_MCD_2001_2005.pdf"

# Extract text from the PDF
pdf_text <- pdf_text(file_path)

# Initialize an empty dataframe to store extracted birth data
birth_data <- data.frame(
  County = character(),
  Municipality = character(),
  Total_Births = integer(),
  Less_1500_Grams = integer(),
  Grams_1500_2499 = integer(),
  Grams_2500_plus = integer(),
  Unknown = integer(),
  stringsAsFactors = FALSE
)

# Define the pattern to detect county names (assuming they appear alone and in uppercase)
county_pattern <- "^[A-Z][A-Za-z\\s]+$"

# Process each page of the PDF
for (page_num in seq_along(pdf_text)) {
  page <- pdf_text[[page_num]]
  lines <- str_split(page, "\n")[[1]]
  county <- NA  # Reset county for each page
  
  for (line_num in seq_along(lines)) {
    line <- str_trim(lines[[line_num]])
    
    # Debugging: Print detected county names
    if (grepl(county_pattern, line) && !grepl("Total", line)) {
      print(paste("Detected County:", line, "on page", page_num))
      county <- line  # Update the current county
      next  # Move to the next line
    }
    
    # Municipality data lines: Assume format "Municipality Name 1,234 12 34 5,678 7"
    if (!is.na(county) && grepl("^\\D+\\s+\\d+", line)) {
      # Remove commas from numbers
      line <- gsub(",", "", line)
      
      parts <- str_match(line, "^(.*?)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)$")
      
      if (!is.na(parts[1, 1])) {
        # Successfully extracted data
        birth_data <- rbind(birth_data, data.frame(
          County = county,
          Municipality = str_trim(parts[1, 2]),
          Total_Births = as.integer(parts[1, 3]),
          Less_1500_Grams = as.integer(parts[1, 4]),
          Grams_1500_2499 = as.integer(parts[1, 5]),
          Grams_2500_plus = as.integer(parts[1, 6]),
          Unknown = as.integer(parts[1, 7]),
          stringsAsFactors = FALSE
        ))
      } else {
        # Debugging: Print lines where extraction failed
        print(paste("Failed to extract municipality data on page", page_num, 
                    "line", line_num, ":", line))
      }
    }
  }
}

# Display the first few rows of the extracted data
head(birth_data)

# Add Wheatland Borough in Mercer to Hermitage City in Mercer
birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] = birth_data[which(birth_data$Municipality == "Hermitage City"),3:7] +
  birth_data[which(birth_data$Municipality == "Wheatland Borough"),3:7]
# remove Wheatland Borough
birth_data = subset(birth_data, Municipality != "Wheatland Borough")
summary(birth_data)

# # Remove any municipalities that are "Total " + word (e.g. "Total Potter")
# birth_data <- birth_data[!grepl("^Total\\s", birth_data$Municipality),]
# write.csv(birth_data, "../data/birth_data.csv", row.names = FALSE)
# 
# # Drop the rows where Total_Births = 1
# birth_data <- read.csv("../data/birth_data.csv")
# #birth_data <- birth_data[birth_data$Total_Births > 1,]

# Calculate the low birth weight rate
birth_data$Low_Birthweight <- birth_data$Less_1500_Grams + birth_data$Grams_1500_2499
birth_data$Num_Births <- birth_data$Total_Births - birth_data$Unknown
birth_data <- select(birth_data, County, Municipality, Low_Birthweight, Num_Births)

# Remove any periods or apostrophes
birth_data$Municipality <- gsub("[[:punct:]]", "", birth_data$Municipality)
# Remoe the space in Dubois
birth_data$Municipality <- gsub("Du Bois City", "DuBois City", birth_data$Municipality)
# Add a space in Bradfordwoods
birth_data$Municipality <- gsub("Bradfordwoods Borough", "Bradford Woods Borough", birth_data$Municipality)
# Remove space in Black Lick
birth_data$Municipality <- gsub("Black Lick Township", "BlackLick Township", birth_data$Municipality)
# Add space to Nantyglo
birth_data$Municipality <- gsub("Nantyglo Borough", "Nanty Glo Borough", birth_data$Municipality)
# Add space to Oilcreek
birth_data$Municipality <- gsub("Oilcreek Township", "Oil Creek Township", birth_data$Municipality)
# Add hyphen to Wilkesbarre
birth_data$Municipality <- gsub("Wilkesbarre City", "Wilkes-Barre City", birth_data$Municipality)
# Add hyphen to Wilkesbarre
birth_data$Municipality <- gsub("Wilkesbarre Township", "Wilkes-Barre Township", birth_data$Municipality)
# Change St Clair Borough to to Saint Clair Borough
birth_data$Municipality <- gsub("St Clair Borough", "Saint Clair Borough", birth_data$Municipality)

# Merge with PA shapefile.
# Source: https://www.pasda.psu.edu/uci/DataSummary.aspx?dataset=41
pamun <- st_read('../shapefiles/PaMunicipalities2025_01.shp')
# to each municipal1 name, add township for CLASS_OF_M = 1TWP 2TWP, borough for BORO, City for CITY, town for TOWN 
pamun$MUNICIPAL1 <- ifelse((pamun$CLASS_OF_M == "1TWP" | pamun$CLASS_OF_M == "2TWP"), 
                           paste0(pamun$MUNICIPAL1, " TOWNSHIP"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "BORO", 
                           paste0(pamun$MUNICIPAL1, " BOROUGH"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "CITY",
                           paste0(pamun$MUNICIPAL1, " CITY"), pamun$MUNICIPAL1)
pamun$MUNICIPAL1 <- ifelse(pamun$CLASS_OF_M == "TOWN",
                           paste0(pamun$MUNICIPAL1, " TOWN"), pamun$MUNICIPAL1)

mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1) # 96.8%
birth_data[which(!toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1),]

# Some of the non-matches are due to "Mount" in the birth_data Muncipality name but "MT" in the shapefile.
# Fix this by replacing "Mount" with "MT" in the birth_data Municipality name, just for these non-matches.
birth_data$Municipality <- ifelse(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1, 
                                  birth_data$Municipality, 
                                  gsub("Mount", "MT", birth_data$Municipality))
birth_data$Municipality[birth_data$Municipality == "Murrysville Municipality"] = "MURRYSVILLE BOROUGH"
birth_data$Municipality[birth_data$Municipality == "Monroeville Municipality" & birth_data$County == 'Allegheny'] = "MONROEVILLE BOROUGH"
birth_data$Municipality[birth_data$Municipality == "Bethel Park Municipality" & birth_data$County == 'Allegheny'] = "BETHEL PARK BOROUGH"

mean(toupper(birth_data$Municipality) %in% pamun$MUNICIPAL1) # 97.1%
birth_data$Municipality <- toupper(birth_data$Municipality)

mean(toupper(birth_data$County) %in% pamun$COUNTY_NAM) # YAY! All match.
birth_data$County <- toupper(birth_data$County)

# Which pamun$MUNICIPAL1 are not in birth_data$Municipality
st_drop_geometry(pamun[which(!(pamun$MUNICIPAL1 %in% birth_data$Municipality)),c('MUNICIPAL1', 'COUNTY_NAM')]) # THERE ARE 14 missing!
# I can't find bw for these.

pamun <- left_join(pamun, birth_data, by = c("COUNTY_NAM" = "County", "MUNICIPAL1" = "Municipality"))

# # Plot low birth weight rate on map using geom_sf
ggplot(pamun) +
  geom_sf(aes(fill = pamun$Low_Birthweight/pamun$Num_Births), color = "black") +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in Pennsylvania")


################################################ Low birth weight in New Jersey ################################################
# Source: https://www-doh.nj.gov/doh-shad/query/builder/birth/BW/LBW.html
bw <- read_xlsx('../data/nj_2001_2005_lowbwt.xlsx') 
# Select Mother's Municipality, Number of Low Birthweight, and Number of Live Births
bw <- bw[,c(2,3,7,8)]
colnames(bw) <- c('Municipality', 'Year', 'Low_Birthweight', 'Num_Births')
# Low Birth weight rate
# Remove commas from numbers
bw$Low_Birthweight <- as.numeric(str_remove_all(bw$Low_Birthweight, ','))
bw$Num_Births <- as.numeric(str_remove_all(bw$Num_Births, ','))
# Convert both to numeric
bw$Low_Birthweight <- as.numeric(bw$Low_Birthweight)
bw$Num_Births <- as.numeric(bw$Num_Births)
# Aggregate data by municipality
bw <- bw %>%
  group_by(Municipality) %>%
  summarise(Low_Birthweight = sum(Low_Birthweight, na.rm = TRUE),
            Num_Births = sum(Num_Births, na.rm = TRUE)) %>%
  ungroup()
# If municipality contains the word "unknown", drop the row
bw <- bw[!grepl('unknown', bw$Municipality, ignore.case = TRUE),]
bw$County <- ifelse(grepl('\\(.*\\)', bw$Municipality),
                    str_extract(bw$Municipality, '\\(.*\\)'),
                    NA)
# Remove parentheses from County
bw$County <- str_remove_all(bw$County, '\\(|\\)')
# Remove the county and the extra space from the municipality name
bw$Municipality <- str_remove_all(bw$Municipality, '\\(.*\\)')
bw$Municipality <- str_trim(bw$Municipality)

# Convert Lyndhurst Borough to Township
bw$Municipality <- ifelse(bw$Municipality == 'Lyndhurst Borough',
                          'Lyndhurst Township',
                          bw$Municipality)
# Convert Caldwell Boro Township to Caldwell Borough
bw$Municipality <- ifelse(bw$Municipality == 'Caldwell Boro Township',
                          'Caldwell Borough',
                          bw$Municipality)
# Convert Essex Fells Township to Essex Fells Borough
bw$Municipality <- ifelse(bw$Municipality == 'Essex Fells Township',
                          'Essex Fells Borough',
                          bw$Municipality)
# Convert Glen Ridge Boro Township to Glen Ridge Borough
bw$Municipality <- ifelse(bw$Municipality == 'Glen Ridge Boro Township',
                          'Glen Ridge Borough',
                          bw$Municipality)
# Convert Orange City to City of Orange Township
bw$Municipality <- ifelse(bw$Municipality == 'Orange City',
                          'City of Orange Township',
                          bw$Municipality)
# Convert Princeton Borough to Princeton
bw$Municipality <- ifelse(bw$Municipality == 'Princeton Borough',
                          'Princeton',
                          bw$Municipality)
# Convert Lavalette Borough to Lavallette Borough
bw$Municipality <- ifelse(bw$Municipality == 'Lavalette Borough',
                          'Lavallette Borough',
                          bw$Municipality)
# Convert Peapack and Gladstone Borough to Peapack-Gladstone Borough
bw$Municipality <- ifelse(bw$Municipality == 'Peapack and Gladstone Borough',
                          'Peapack-Gladstone Borough',
                          bw$Municipality)
# Convert Avon-By-The-Sea Borough to Avon-by-the-Sea Borough
bw$Municipality <- ifelse(bw$Municipality == 'Avon-By-The-Sea Borough',
                          'Avon-by-the-Sea Borough',
                          bw$Municipality)

# Read in NJ Municipality shapefiles
# https://njogis-newjersey.opendata.arcgis.com/datasets/municipal-boundaries-of-nj-hosted-3424/explore
nj_mun <- st_read('../shapefiles/NJ_Municipal_Boundaries/NJ_Municipal_Boundaries_3424.shp')
mean(bw$Municipality %in% nj_mun$MUN_LABEL) # 0.986
mean(nj_mun$MUN_LABEL %in% bw$Municipality) # 1

#bw$Low_Birthweight_Rate <- bw$Low_Birthweight / bw$Num_Births
#hist(bw$Low_Birthweight_Rate)

nj_mun$Low_Birthweight <- NA
nj_mun$Num_Births <- NA
for (i in 1:nrow(bw)){
  if (bw$Municipality[i] %in% nj_mun$MUN_LABEL){
    # If one match for municipality name, 
    if (is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$Municipality[i])
      nj_mun$Low_Birthweight[j] <- bw$Low_Birthweight[i]
      nj_mun$Num_Births[j] <- bw$Num_Births[i]
    }
    # If two matches for municipality name, choose the one with the same county
    if (!is.na(bw$County[i])){
      j <- which(nj_mun$MUN_LABEL == bw$Municipality[i] & nj_mun$COUNTY == toupper(bw$County[i]))
      nj_mun$Low_Birthweight[j] <- bw$Low_Birthweight[i]
      nj_mun$Num_Births[j] <- bw$Num_Births[i]
    }
  }
}

# Plot low birth weight rate on map using geom_sf
ggplot(nj_mun) +
  geom_sf(aes(fill = nj_mun$Low_Birthweight/nj_mun$Num_Births), color = NA) +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in New Jersey")

# Combine the PA shapefile with the NJ shapefile
pamun <- st_transform(pamun, crs = st_crs(nj_mun))
pamun <- select(pamun, COUNTY, MUNICIPAL1, Low_Birthweight, Num_Births)
nj_mun <- select(nj_mun, COUNTY, MUN_LABEL, Low_Birthweight, Num_Births)
# Rename MUNICIPAL1 and MUN_LABEL each to Municipality
colnames(pamun)[2] <- 'Municipality'
colnames(nj_mun)[2] <- 'Municipality'
# Combine the two dataframes
mun_baseline <- rbind(pamun, nj_mun)

# Plot low birth weight rate on map using geom_sf
ggplot(mun_baseline) +
  geom_sf(aes(fill = mun_baseline$Low_Birthweight/mun_baseline$Num_Births), color = NA) +
  scale_fill_viridis_c(name = "Low Birthweight Rate") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Low Birthweight Rate by Municipality in New Jersey and Pennsylvania")
st_write(mun_baseline, '../shapefiles/Municipality_Low_Birthweight_Rate_2001_2005.shp')

