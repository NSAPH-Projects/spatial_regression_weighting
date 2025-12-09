# For chemicals and contaminants
library(readxl)
library(dplyr)
library(PeriodicTable)

# https://www.epa.gov/superfund/superfund-data-and-reports
# Read in contaminants of concern xlsx

url <- 'https://semspub.epa.gov/src/document/HQ/406203'
tmp <- tempfile(fileext = ".xlsx")

download.file(url, destfile = tmp, mode = "wb")
chemical_medium <- read_xlsx(tmp)

# Remove first row, assign second row to column names
colnames(chemical_medium) <- as.character(chemical_medium[1, ])

chemical_medium <- chemical_medium %>%
  mutate(`Contaminant Name` = toupper(`Contaminant Name`))

sort(table(chemical_medium$`Contaminant Name`))
sort(table(chemical_medium$`Media`))

medium_groups <- list(
  Gas = c("Air", "Landfill Gas", "Soil Gas"),
  Solid = c("Soil", "Sediment", "Sludge", "Solid Waste",
                 "Residuals", "Debris", "Buildings/Structures"),
  Water = c("Groundwater", "Surface Water", "Leachate",
            "Liquid Waste", "Free-phase NAPL"),
  Other = c("Other", "Fish Tissue")
)

# Get the elements from PeriodTable which are metals
metals <- periodicTable$name[periodicTable$type %in% c('Alkali Metal', 'Alkaline Earth Metal',
                                                        'Transition Metal', 'Transactinide',
                                                        'Lanthanide', 'Actinide', 'Metal')]
is_metal <- chemical_medium$`Contaminant Name` %in% toupper(metals)

# Define vector of chlorinated solvents and VOCs (based on EPA/ATSDR examples)
chlorinated_vocs <- c(
  "TRICHLOROETHENE",
  "TETRACHLOROETHENE",
  "1,1,1-TRICHLOROETHANE",
  "1,1,2-TRICHLOROETHANE",
  "TRICHLOROETHANE (MIXED ISOMERS)",
  "1,1-DICHLOROETHENE",
  "1,2-DICHLOROETHENE (CIS AND TRANS MIXTURE)",
  "CIS-1,2-DICHLOROETHENE",
  "TRANS-1,2-DICHLOROETHENE",
  "1,2-DICHLOROETHANE",
  "CHLOROFORM",
  "CARBON TETRACHLORIDE",
  "DICHLOROMETHANE (METHYLENE CHLORIDE)",
  "METHYLENE CHLORIDE",
  "CHLOROETHENE (VINYL CHLORIDE)",
  "VINYL CHLORIDE",
  "CHLOROBENZENE",
  "1,2-DICHLOROBENZENE",
  "1,3-DICHLOROBENZENE",
  "1,4-DICHLOROBENZENE",
  "CHLOROETHANE",
  "1,2-DICHLOROPROPANE",
  "1,2,3-TRICHLOROPROPANE",
  "PENTACHLOROETHANE",
  "1,1,2,2-TETRACHLOROETHANE",
  "1,1,1,2-TETRACHLOROETHANE",
  "TRICHLOROFLUOROMETHANE",
  "DICHLORODIFLUOROMETHANE",
  "CHLOROMETHANE",
  "BROMODICHLOROMETHANE",
  "DIBROMOCHLOROMETHANE",
  "BROMOFORM"
)
is_voc <- chemical_medium$`Contaminant Name` %in% chlorinated_vocs

# Persistent Organic Pollutants
# According to Stockholm Convention https://www.pops.int/TheConvention/ThePOPs/AllPOPs/tabid/2509/Default.aspx
# And renaming given what's in my data
pops <- c(
  "ALDRIN",
  "ALPHA-HEXACHLOROCYCLOHEXANE",
  "BETA-HEXACHLOROCYCLOHEXANE",
  "CHLORDANE",
  "CHLORDECONE",
  "DDT",
  "DICOFOL",
  "DIELDRIN",
  "ENDOSULFAN",
  "ENDRIN",
  "HEPTACHLOR",
  "HEXABROMOBIPHENYL",
  "HEXABROMOCYCLODODECANE",
  "HEXABROMODIPHENYL ETHER",
  "HEPTABROMODIPHENYL ETHER",
  "HEXACHLOROBENZENE",
  "HEXACHLOROBUTADIENE",
  "LINDANE",
  "METHOXYCHLOR",
  "MIREX",
  "PENTACHLOROBENZENE",
  "PENTACHLOROPHENOL",
  "PERFLUOROOCTANE SULFONIC ACID",
  "PERFLUOROOCTANOIC ACID",
  "PERFLUOROHEXANE SULFONIC ACID",
  "POLYCHLORINATED BIPHENYLS (PCBS)",
  "POLYCHLORINATED NAPHTHALENES",
  "SHORT-CHAIN CHLORINATED PARAFFINS",
  "MEDIUM-CHAIN CHLORINATED PARAFFINS",
  "TETRABROMODIPHENYL ETHER",
  "PENTABROMODIPHENYL ETHER",
  "DECABROMODIPHENYL ETHER",
  "DECHLORANE PLUS",
  "TOXAPHENE",
  "UV-328"
)

is_pop <- chemical_medium$`Contaminant Name` %in% pops
chem_mm <- chemical_medium %>%
  mutate(
    Gas   = Media %in% medium_groups$Gas,
    Solid = Media %in% medium_groups$Solid,
    Water = Media %in% medium_groups$Water,
    Metal_flag = is_metal,
    VOC_flag   = is_voc
  )

# Polycyclic Aromatic Hydrocarbons (PAHs)
# https://www.ncbi.nlm.nih.gov/books/NBK217760/
# And what's in my data
pahs <- c(
  "ACENAPHTHENE",
  "ACENAPHTHYLENE",
  "ANTHRACENE",
  "BENZO[A]ANTHRACENE",
  "BENZO[A]PYRENE",
  "BENZO[B]FLUORANTHENE",
  "BENZO[J]FLUORANTHENE",
  "BENZO[K]FLUORANTHENE",
  "BENZO[G,H,I]PERYLENE",
  "CHRYSLENE",
  "DIBENZO[A,H]ANTHRACENE",
  "FLUORANTHENE",
  "FLUORENE",
  "INDENO[1,2,3-CD]PYRENE",
  "NAPHTHALENE",
  "PHENANTHRENE",
  "PYRENE"
)
is_pah <- chemical_medium$`Contaminant Name` %in% pahs

site_flags <- chem_mm %>%
  group_by(`EPA ID`) %>%                     
  summarize(
    metal = as.integer(any(replace_na(Metal_flag, FALSE))),
    voc   = as.integer(any(replace_na(VOC_flag,   FALSE))),
    pop   = as.integer(any(replace_na(is_pop,     FALSE))),
    pah   = as.integer(any(replace_na(is_pah,     FALSE))),
    gas   = as.integer(any(replace_na(Gas,        FALSE))),
    solid = as.integer(any(replace_na(Solid,      FALSE))),
    water = as.integer(any(replace_na(Water,      FALSE))),
    .groups = "drop"
  )
# NOTE: ALL SITES HAD AT LEAST ONE POP OR PAH LISTED! These will have to be ignored as covariates.
 
colnames(site_flags)[1] <- 'S_EPA_I'
write.csv(site_flags, "data/site_chemical_medium_flags.csv", row.names = FALSE)
