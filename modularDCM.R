rm(list=ls())

# DONE
source("/home/flo/Documents/DCM_BlackSea/data_extraction_cleaningV2.R")
source("/home/flo/Documents/DCM_BlackSea/data_fitV2.R") 
source("/home/flo/Documents/DCM_BlackSea/data_densityV2.R")
source("/home/flo/Documents/DCM_BlackSea/data_fit_BBP.R")

# ADD SOME DATA FOR THE .Rmd file
# get unique positions of CHLA profiles
lat_lon <- subset(profiles, select = c("lat","lon","id"))
lat_lon <- lat_lon[!duplicated(lat_lon),]

# TBD
#load("light_iso_argo.Rdata")


# TO REMOVE?
#source("/home/flo/Documents/DCM_BlackSea/DCM_BlackSea/light.R")
#source("/home/flo/Documents/DCM_BlackSea/DCM_BlackSea/test_isolumesV3.R")

save.image(file='data.RData') # save workspace directly

