# Deep Chlorophyll Maximum dynamics in the Black Sea

### Method

- Extract data (CHLA, CDOM, BBP700, DOXY, PAR) and clean them with the associated flags (see ARGO DOCUMENTATION)
- Correct CHLA data based on the paper of Xing et al., 2017 (https://doi.org/10.1002/lom3.10144)
  - Data_extraction_cleaningV2.R
  
- Fit the data with Gaussian curves if possible to spot DCM
  - data_fitV2.R and data-fit_BBP.R for BBP profiles
  
- Compute the density (to investigate the DCM density-related hypothesis (https://doi.org/10.1002/gbc.20093)
  - data_densityV2.R
  
### Small VISU app related to this

See at https://fricour.shinyapps.io/blackseadcm/ (source code embedded in it but it can still be found in this directory with the R markdown file (actually a flex dashboard) ===> XX.Rmd

  
