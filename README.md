# Understanding Spatial Regression Models from a Weighting Perspective in an Observational Study of the Impact of Superfund Cleanups on Birth Weight

This repository contains the code for the simulation and data application sections of the paper. The workflow is:

- `funcs.R` contains all utility functions for the simulation and data application. 
- `simulation` contains the code for the simulation section of the paper. The simulation is fully reproducible. 
- `application` contains the code for the data application section of the paper. 

### Data
We link publicly available U.S. Environmental Protection Agency information on 1,429 Superfund sites with high-resolution 1990 Census covariates and Medicaid birth outcomes from 2016–2018. We define a binary exposure indicating whether a site was remediated between 1991 and 2015. For each site and its 2-km buffer, we assemble eleven sociodemographic covariates, five indicators of contaminant types, and site scores. Exposure and confounders are publicly available on Harvard Dataverse at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EKSCCU. Birth outcomes (the percentages of small vulnerable newborns and congenital anomalies) are aggregated from zip code tabulation area-level Medicaid claims, yielding a final analytic sample of 1,079 Superfund sites. Medicaid claims data is stored on the Harvard Regulated Data Environment and is not publicly available. 

### Data Citations
  - U.S. Census Bureau (1991), 1990 Decennial Census, Summary Tape Files 1 and 3A (STF1 and 3A).
  - Jonathan Schroeder, David Van Riper, Steven Manson, Katherine Knowles, Tracy Kugler, Finn Roberts, and Steven Ruggles. 
        IPUMS National Historical Geographic Information System: Version 20.0 
        [dataset]. Minneapolis, MN: IPUMS. 2025. 
        http://doi.org/10.18128/D050.V20.0
  - U.S. Environmental Protection Agency (2025). Superfund National Priorities List (NPL) Where you Live Map. ArcGIS Feature Service. https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=33cebcdfdd1b4c3a8b51d416956c41f1
  - U.S. Environmental Protection Agency (2025). Contaminant of Concern Data for Decision Documents by Media, FYs 1981-2024 (Final NPL, Deleted NPL, and Superfund Alternative Approach Sites), https://www.epa.gov/superfund/superfund-data-and-reports.
  - U.S. Environmental Protection Agency, Office of Mission Support (2025). NPL Superfund Site Boundaries (EPA Public), https://www.arcgis.com/home/item.html?id=d6e1591d9a424f1fa6d95a02095a06d7.
  - Centers for Medicare & Medicaid Services (2020). Medicaid and CHIP T-MSIS Analytic Inpatient file. **Restricted data containing protected health information; not publicly accessible.**
  
### Acknowledgements
We are grateful to Shreya Nalluri, Amruta Nori-Sarma, Mary Willis, and Flannery Black-Ingersoll for their thoughtful discussions on Medicaid birth outcome data. The simulations in this paper were run on the FASRC Cannon cluster, supported by the FAS Division of Science Research Computing Group at Harvard University. The outcome analysis was run on Regulated Data (ReD) Environment, operated by Harvard University Information Technology and supported by the Office of the Vice Provost for Research at Harvard University. 
 
### Paper
The paper associated with this repository is available on arXiv:
["Understanding Spatial Regression Models from a Weighting Perspective in an Observational Study of Superfund Remediation"](https://arxiv.org/pdf/2508.19572) (Woodward, Dominici, Zubizarreta 2025).

### Workflow for reproducibility
To reproduce the results in the paper, please follow these steps. 

0. Download publicly available data:
   - Download 1990 census tract shapefiles and the following variables from the 1990 Decennial Census via [NHGIS](https://catalog.data.gov/dataset/census-tracts-in-1990/resource/18c7c256-0b30-41e0-887e-5b0b2ce2dc20): total population, population of individuals identifying as Hispanic, Black, Indigenous, or Asian, median household income, median house value, percentage below the poverty line, tenure status, population with a high school diploma, and the median year of housing construction. See codebooks in `data/nhgis0002_csv` for specific variable names. Save all files in the `data/` folder.
1. Run `application/01_chemical_medium_scrape.R` to process chemical and medium information from the contaminants of concern file.
2. Run `application/02_preprocessing.R` to web scrape remediation dates and site scores from EPA Superfund site profiles, extract Superfund site boundary information from EPA ArcGIS, and merge with 1990 Census tract covariates and the chemical/medium information from step 1. After area weighting, this script produces the Superfund site-level dataset of binary exposure and covariates called `buffers` which can be accessed by executing `load('data/preprocessed_superfunds.RData')`.
3. Run `application/03_diagnostics.R` to produce Figures 1, 2, 7, and 8.
4. Run `application/04_balanace_gamma_diagnostics.R` to produce Figure 6. 
5. Run `application/05_random_fixed_effects.R` to produce Figure 5. 
6. Run `application/06_outcome_analysis.R` to produce Figures 3, 4, and 12 as well as Tables 1, 3, 5, 6, and 7.
7. Run `simulation/01_simulation_preprocessing.R` to preprocess the data needed for the simulation study and generate Figure 10.
8. Execute `sbatch submit_jobs.sh` to run the simulation analyses on a computing cluster. This will generate results as CSV files in the `results_Oct30/` folder.
9. Run `simulation/03_analysis.R` to produce Table 2 and Figure 9.

### Contact us
If you have any questions, please contact us at [swoodward@g.harvard.edu](mailto:swoodward@g.harvard.edu).