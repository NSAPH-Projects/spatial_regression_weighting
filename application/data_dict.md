# README: Preprocessed Superfunds Data

## Overview
This dataset contains preprocessed treatment and confounder data for ["Understanding Spatial Regression Models from a Weighting
Perspective in an Observational Study of Superfund Remediation"](https://arxiv.org/pdf/2508.19572) (Woodward, Dominici, Zubizarreta 2025).  
The file `preprocessed_superfunds.RData` includes exposure, confounder, distance, and adjacency information for 1,583 Superfund sites listed on the National Priorities List in 2001.

---
            
## Contents of `preprocessed_superfunds.RData`

### 1. `X` (Design matrix)
A data frame of covariates for each buffer around a Superfund site.

| Variable | Type | Units | Description |
|----------|------|-------|-------------|
| Intercept | numeric | – | Constant 1. |
| Site_Score | numeric | – | NPL site score from EPA. |
| metal | integer (0/1) | – | Indicator for presence of metals among contaminants of concern at the site. |
| voc | integer (0/1) | – | Indicator for presence of volatile organic compounds among contaminants of concern at the site. |
| gas | integer (0/1) | – | Indicator for contamination through air/soil/landfill gas at the site. |
| solid | integer (0/1) | – | Indicator for contamination through solid media (soil, sediment, sludge, solid waste, etc) at the site. |
| water | integer (0/1) | – | Indicator for contamination through liquid media (groundwater, surface water, leachate, etc) at the site. |
| population_density | numeric | people per square mile | Estimated population density within 2-km buffer (from Census 1990 tracts, area-weighted). |
| percent_hispanic | numeric | proportion (0–1) | Area-weighted mean Hispanic proportion. |
| percent_black | numeric | proportion (0–1) | Area-weighted mean Black proportion. |
| percent_indigenous | numeric | proportion (0–1) | Area-weighted mean American Indian/Alaska Native proportion. |
| percent_asian | numeric | proportion (0–1) | Area-weighted mean Asian proportion. |
| percent_renter_occupied | numeric | proportion (0–1) | Area-weighted mean proportion of renter-occupied housing units. |
| median_household_income | numeric | US dollars | Area-weighted median household income. |
| median_house_value | numeric | US dollars | Area-weighted median house value. |
| percent_poverty | numeric | proportion (0–1) | Area-weighted mean poverty rate. |
| percent_high_school_grad | numeric | proportion (0–1) | Area-weighted mean of tract education levels (HS graduate or higher). |
| median_year_structure_built | numeric | year | Area-weighted median housing unit year built. |

---

### 2. `Z` (Treatment indicator)
- **Type**: integer (0/1) vector  
- **Definition**: `1` if the site was cleaned up (deleted from NPL) between 1991–2015, `0` otherwise.  
- **Source**: EPA NPL status table.

---

### 3. `dmat` (Distance matrix)
- **Type**: numeric matrix (n × n), symmetric with zeros on diagonal  
- **Units**: kilometers  
- **Definition**: Pairwise distances between site centroids (Haversine).  

---

### 4. `clusters`
- **Type**: integer vector  
- **Definition**: State-level numeric encoding.  

---

### 5. `adjacency_matrix`
- **Type**: sparse symmetric 0/1 `Matrix`  
- **Definition**: Five-nearest-neighbors graph; `(i,j)=1` if `j` is among the five nearest sites to `i` (symmetrized).  

---

### 6. `buffers` (sf object)
An `sf` and `data.frame` object of 2-km buffer polygons around Superfund site areas with merged covariates.

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| Site_EPA_ID | character | – | EPA site identifier. |
| Site_Score | numeric | – | NPL site score. |
| metal | integer (0/1) | – | Indicator for presence of metals among contaminants of concern at the site. |
| voc | integer (0/1) | – | Indicator for presence of volatile organic compounds among contaminants of concern at the site. |
| gas | integer (0/1) | – | Indicator for presence of contamination through air/soil/landfill gas at the site. |
| solid | integer (0/1) | – | Indicator for presence of contamination through solid media (soil, sediment, sludge, solid waste, etc) at the site. |
| water | integer (0/1) | – | | Indicator for presence of contamination through liquid media (groundwater, surface water, leachate, etc) at the site. |
| Z | integer (0/1) | – | Treatment indicator (see above). |
| cluster | integer | – | Cluster ID (see above). |
| Latitude | numeric | degrees | Site latitude. |
| Longitude | numeric | degrees | Site longitude. |
| population_density | numeric | people per square mile | Estimated population density within 2-km buffer (from Census 1990 tracts, area-weighted). |
| percent_hispanic | numeric | proportion (0–1) | Area-weighted mean Hispanic proportion. |
| percent_black | numeric | proportion (0–1) | Area-weighted mean Black proportion. |
| percent_white | numeric | proportion (0–1) | Area-weighted mean White proportion. |
| percent_indigenous | numeric | proportion (0–1) | Area-weighted mean American Indian/Alaska Native proportion. |
| percent_asian | numeric | proportion (0–1) | Area-weighted mean Asian proportion. |
| percent_renter_occupied | numeric | proportion (0–1) | Area-weighted mean proportion of renter-occupied housing units. |
| median_household_income | numeric | US dollars | Area-weighted median household income. |
| median_house_value | numeric | US dollars | Area-weighted median house value. |
| percent_poverty | numeric | proportion (0–1) | Area-weighted mean poverty rate. |
| percent_high_school_grad | numeric | proportion (0–1) | Area-weighted mean of tract education levels (HS graduate or higher). |
| median_year_structure_built | numeric | year | Area-weighted median housing unit year built. |
| geoms | sfc_POLYGON | meters | 2-km buffer polygon geometry (EPSG:5070). |

---

## Data Sources and Processing
- **EPA NPL data**: U.S. Environmental Protection Agency, National Priorities List (NPL) sites.  
- **1990 Decennial Census data**: U.S. Census Bureau (SF1 and SF3 tables) accessed via NHGIS.  
- **Processing steps**:
  - Constructed 2-km buffers around site areas.  
  - Merged chemical/medium indicators from contaminants of concern.
  - Intersected buffers with Census tracts (1990).  
  - Derived area-weighted averages for proportions and area-weighted medians for medians.  
  - Computed population densities from tract-level total population and area.  
  - Built distance and adjacency matrices for spatial modeling.  

---

Final sample size: 1,429 Superfund sites.