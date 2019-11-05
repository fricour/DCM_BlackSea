#############
# LIBRARIES #
#############

library(ncdf4)
library(ggplot2)
library(plyr)
library(reshape2)
library(chron)
library(viridis)
library(TTR)
library(minpack.lm)
library(nls2)
library(nlstools)
library(data.table)
library(gsw)
library(gridExtra)
library(car)
library(oce)
library(ggmap)
library(ggalt)
library(RColorBrewer)
library(lubridate)
library(lattice)
library(suncalc)
library(dplyr)
library(plyr)
library(ggplot2)

#############
# FUNCTIONS #
#############

#1 Extract data from any variable available in the float data files
ExtractVar <- function(Var,FloatInfo){
  
  with(FloatInfo,{
    # This function should return a dataframe for the chosen variable with value, qc, iprofile and ilevel
    lvar             <- ncvar_get(ncfile,Var)
    lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
    lvar_qctab       <- llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab <- do.call(cbind,lvar_qctab)
    
    lvar_dir       <- ncvar_get(ncfile,"DIRECTION")
    lvar_direction <- llply(lvar_dir,function(dirstring){
      strsplit(dirstring,split="")
    })
    lvar_direction <- unlist(lvar_direction)
    # making dataframes, removing the NANs  
    alevels <- 1:N_LEVELS
    d <- ldply(as.list(1:N_PROF),function(iprof){
      indexes <- !(is.na(lvar[,iprof])|is.na(lvar[,iprof]))
      if(sum(indexes) == 0){
        return (data.frame())
      }
      
      data.frame(value    = lvar[indexes,iprof],
                 qc       = as.integer(lvar_qctab[indexes,iprof]),
                 alevel   = alevels[indexes],
                 depth    = pres[indexes,iprof],
                 dir      = lvar_direction[iprof],
                 aprofile = iprof,
                 variable = Var)
    })
    
    d$juld <- juld[d$aprofile]
    d$lon  <- lon[d$aprofile]
    d$lat  <- lat[d$aprofile]
    
    return(d=d)
  })
}  

#2 five-points median filter
mmed <- function(x,n=5){runmed(x,n)}

#3 Determination of day/night profiles
DayOrNight <- function(juld, lat, lon){
  origin <- ymd_hms("1950-01-01 00:00:00") # julian day origin for the Argo program
  juld <- origin + juld * 3600 * 24 #conversion of juld in working day and hour
  #Below we split the julian day into day - month - year - hour - minute
  day <- format(as.POSIXct(strftime(juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%d")
  month <- format(as.POSIXct(strftime(juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%m")
  year <- format(as.POSIXct(strftime(juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%Y")
  hour <- format(as.POSIXct(strftime(juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%H")
  minute <- format(as.POSIXct(strftime(juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%M") 
  date <- as.Date(as.POSIXct(juld, 'Europe/Sofia')) # convert juld into a a "Date" object
  data <- data.frame(date = date, lat = lat, lon = lon) # create a dataframe for the function getSunlightTimes
  sun_cycle <- getSunlightTimes(data = data, keep = c("sunset", "sunrise"), tz = "Europe/Sofia")# get the sunset and sunrise hours
  
  # local sunrise/sunset hours
  hour_sunrise <- as.numeric(format(sun_cycle$sunrise, "%H"))
  hour_sunset <- as.numeric(format(sun_cycle$sunset, "%H"))
  
  if(as.numeric(hour) < hour_sunset & as.numeric(hour) > hour_sunrise){
    sun <- "day"
  }else{
    sun <- "night"
  }
  return(sun)
}

#4 Function to reinitialize dataframes after cleaning (row numbers, new IDs)
clean_remove <- function(dataframe, vector_filter){
  dataframe <- dataframe[dataframe$juld %in% vector_filter,]
  rownames(dataframe) <- NULL
  dataframe <- transform(dataframe, id=as.numeric(factor(juld)))
  dataframe <- dataframe[order(dataframe$id),]
}

#5 Function to normalize data
normalize <- function(data){
  data <- (data - min(data, na.rm = T))/(max(data, na.rm = T) - min(data, na.rm = T))
}

##########
# FLOATS #
##########

filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")
#filename <- c("6900807_Sprof.nc","6901866_Sprof.nc","7900591_Sprof.nc","7900592_Sprof.nc") #let's use the synthetic profiles?

#####################
# GET PROFILES DATA #
#####################

# Acquisition of CHLA and BBP (non adjusted !) and CDOM (if available) profiles from each float
profiles <- ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = F, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  
  FloatInfo <- list(N_PROF=N_PROF,
                    N_LEVELS=N_LEVELS,
                    juld=juld,
                    pres=pres,
                    lon=lon,
                    lat=lat)
  
  chladf <- ExtractVar("CHLA", FloatInfo)
  bbpdf <- ExtractVar("BBP700", FloatInfo)
  
  if (file == "7900591_Mprof.nc" | file == "7900592_Mprof.nc"){ #change M for merged or S for synthetic profiles
    cdomdf <- as.data.frame(rep(NA, length(chladf$value)))
    colnames(cdomdf) <- "value"
  }else{
    cdomdf <- ExtractVar("CDOM", FloatInfo)
  }
  
  # for a dimension mismatch issue occurring with FLOAT WMO 6900807 (between BBP and CHLA)
  if(file == "6900807_Mprof.nc"){
    chladf <- join(bbpdf, chladf, by = c("alevel", "aprofile"))
    chladf <- chladf[,11:18]
    cdomdf <- join(bbpdf, cdomdf, by = c("alevel", "aprofile"))
    cdomdf <- cdomdf[,11:18]
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  day_night <- DayOrNight(chladf$juld[1], chladf$lat[1], chladf$lon[1])
  
  data.frame(depth    = chladf$depth,
             juld     = chladf$juld,
             fluo     = chladf$value,# fluo = chla converted from fluo data
             bbp      = bbpdf$value,
             cdom     = cdomdf$value,
             qc_fluo  = chladf$qc,
             qc_bbp   = bbpdf$qc, #no qc_cdom because qc_cdom is always 0 so we cannot exploit it
             day      = month.day.year(chladf$juld,c(1,1,1950))$day,
             month    = month.day.year(chladf$juld,c(1,1,1950))$month,
             year     = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = chladf$lon,
             lat      = chladf$lat,
             dir      = chladf$dir,
             day_night = day_night,
             platform = as.numeric(unique(id))
  )
})

######################
### DATA CLEANING  ###
######################

# REMOVE DATA OF 2013 AND 2019 (incomplete years)
profiles <- profiles[!(profiles$year == 2013 | profiles$year == 2019),]
# REMOVE BAD DATA (I.E. QC = 4 FOR NEGATIVE SPIKES, JUMPS, ETC.)
profiles <- profiles[-(which(profiles$qc_fluo == 4)),] # bad chla data
profiles <- profiles[-(which(profiles$qc_bbp %in% c(3,4))),] # bad bbp data
# REMOVE DESCENT PROFILES
profiles <- profiles[-(which(profiles$dir == "D")),]
# REMOVE PROFILES WITH NO LAT/LON
profiles <- profiles[-(which(is.na(profiles$lat == T))),]

# CREATION OF PROFILE IDs & REORDER DATA FRAME ACCORDING TO IT
profiles <- clean_remove(profiles, unique(profiles$juld))

# MEDIAN FILTER TO REMOVE POSITIVE SPIKES (THEY MAY RETAIN INFORMATION BUT STRONG SPIKES WILL NOT BE RETAINED FOR THIS STUDY)
smoothed_fluo <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

smoothed_bbp <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  tmp <- mmed(tmp$bbp, 5)
  data.frame(smoothed_bbp = tmp)
})

smoothed_cdom <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  if (tmp$platform[1] == 6900807 | tmp$platform[1] == 6901866){
    tmp <- mmed(tmp$cdom, 5)
  }else{
    tmp <- tmp$cdom
  }
  data.frame(smoothed_cdom = tmp)
})

# REPLACE FLUO WITH SMOOTHED FLUO
profiles[,3] <- smoothed_fluo
profiles[,4] <- smoothed_bbp
profiles[,5] <- smoothed_cdom

# remove data not needed anymore (free some space)
rm(smoothed_bbp, smoothed_fluo, smoothed_cdom, ncfile)

# Get some basic info on each profile
argodf <- ddply(profiles,~juld,summarize,
               depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],# to not take into account surface depth where chlorophyll can be minimum
               bottomdepth = max(depth),
               min_depth = min(depth),
               platform = platform[1],
               month = month[1],
               year = year[1], DOY = DOY[1],
               n = length(fluo)) #n: number of points in the profile

# IDs
argodf <- clean_remove(argodf, unique(argodf$juld))

# Get profiles id that do not answer to some criteria
crit1 <- argodf[(argodf$min_depth > 5),]# all stuck profiles except one (add QC on pressure to get rid of them directly)
crit2 <- argodf[(argodf$bottomdepth < 40),] # 50m works as well, justification: le dcm est en-dessous de 50m 

# remove a profile (.6 24895.96 0.8784 0.00916 2.7240 2 2 2 2018 59 34.85083 44.71711 A night 6900807 652) that shows a DCM at > 150m deep (cannot correct with FDOM method) -> remove it
crit3 <- argodf[argodf$platform == 6900807 & argodf$DOY == 59 & argodf$year == 2018,]

criteriadf <- rbind(crit1, crit2, crit3)
criteriadf <- criteriadf[(!duplicated(criteriadf)),] # this line can be removed because it does not change anything for our case

#remove profiles that have been categorized as bad profiles
argodf <- argodf[-(which(argodf$id %in% criteriadf$id)),]
profiles <- profiles[-which(profiles$id %in% criteriadf$id),]
profiles <- clean_remove(profiles, unique(profiles$juld))

# clean workspace
rm(crit1, crit2, crit3, criteriadf)
  
###################################
# COMPUTATION OF DENSITY PROFILES #
###################################

density_profiles <- ldply(as.list(filename),function(file){
  
  ncfile   <<- nc_open(file, write = FALSE, verbose = F, suppress_dimvals = FALSE)
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  FloatInfo<-list(N_PROF=N_PROF,
                  N_LEVELS=N_LEVELS,
                  juld = juld,
                  pres=pres,
                  lon=lon,
                  lat=lat)
  
  #NO ADJUSTED VALUE FOR PSAL AND TEMP FOR EACH FLOAT
  tempdf <- ExtractVar("TEMP",FloatInfo)
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- ExtractVar("PSAL",FloatInfo)
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA
  
  #DENSITY ANOMALY BASED ON TEOS-10
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,tempdf$depth)
  sigma <- gsw_sigma0(psal,temp)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  data.frame(sigma = sigma,
             depth = tempdf$depth, 
             juld  = tempdf$juld, 
             dir   = tempdf$dir,
             lon   = tempdf$lon,
             lat   = tempdf$lat,
             #day      = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month    = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year     = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY   = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                         format ="%j")),
             platform = as.numeric(unique(id)))
})

# New IDs
rownames(profiles) <- NULL
density_profiles <- clean_remove(density_profiles, unique(profiles$juld))
argodf <- clean_remove(argodf, unique(profiles$juld))

#COMPUTE MLD (criteria can thus be modified according to the definition of MLD that you choose)
sigma_criteria <- 0.03
depth_ref <- 10

MLDdf <- ldply(as.list(1:length(unique(density_profiles$id))), function(i){
  
  tmp <- density_profiles[density_profiles$id == i,]
  rownames(tmp) <- NULL
  
  if(all(is.na(tmp$sigma)) == TRUE){
    MLD <- NA
    sigma_surface <- NA
  }else if(length(tmp$sigma[!is.na(tmp$sigma)==TRUE])>=2){
    sigma_surface<-NA
    sigma_surface <- approx(tmp$depth,tmp$sigma,depth_ref)$y
  }
  
  if(is.na(sigma_surface) == FALSE){
    MLD <- max(tmp$depth[tmp$sigma <= (sigma_surface + sigma_criteria)], na.rm = T)
  }else{
    MLD <- NA
  }
  
  if(is.na(MLD) == TRUE){
    sigmaMaxMLD <- NA
  }else{
    sigmaMaxMLD <- mean(tmp$sigma[which(tmp$depth == MLD)], na.rm =T)
  }
  
  show(i)
  
  data.frame(sigma_surface = sigma_surface, MLD = MLD, sigmaMaxMLD = sigmaMaxMLD,
             juld = tmp$juld[1], lat=tmp$lat[1], lon=tmp$lon[1], day = month.day.year(tmp$juld[1],c(1,1,1950))$day,
             month = month.day.year(tmp$juld[1],c(1,1,1950))$month,
             year = month.day.year(tmp$juld[1],c(1,1,1950))$year,
             DOY = as.integer(strftime(as.Date(tmp$juld[1],origin = '1950-01-01'), 
                                       format ="%j")),
             platform = tmp$platform[1])
})

#REMOVE BAD PROFILES (NO MLD DATA AVAILABLE)
bad_mld_juld <- MLDdf$juld[which(is.na(MLDdf$MLD) == TRUE)]
profiles <- profiles[-which(profiles$juld %in% bad_mld_juld),]
profiles <- clean_remove(profiles, unique(profiles$juld))

# reorder dataframes
density_profiles <- clean_remove(density_profiles, unique(profiles$juld))
argodf <- clean_remove(argodf, unique(profiles$juld))
MLDdf <- clean_remove(MLDdf, unique(profiles$juld))

# add a column to the 'profiles' dataframe with OLD_CHLA (will save the uncorrected chla, can be useful for some plots)
profiles$old_fluo <- profiles$fluo

## In reality, it could be a script based not on the ID but on the JULD (avoid all the rownames init and so forth.. To meditate)

#FDOM-BASED CORRECTION AND MINIMUM-OFFSET CORRECTION
FDOM_OR_MINIMUM_OFFSET_CORRECTION <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id == i,]
  MaxDepth <- argodf$bottomdepth[i]#Max depth of the profile
  TopDepth <- argodf$depthmin[i]#Apparent minimum 
  
  if(all(is.na(tmp$cdom)) == TRUE){#no cdom sensor hence minimum-offset correction procedure
    offset <- tmp$fluo[which(tmp$depth == TopDepth)]
    fluo_cor <- tmp$fluo - offset
    fluo_cor[which(tmp$depth == TopDepth):length(tmp$depth)] <- 0
    slope_fdom <- NA
    C <- NA
  }else{#FDOM-based method
    calibrange <- tmp[which(tmp$depth == TopDepth):length(tmp$depth),]
    linearMod <- lm(fluo ~ cdom, data = calibrange)
    slope_fdom <- coef(linearMod)[[2]]
    C <- coef(linearMod)[[1]]
    fluo_cor <- tmp$fluo - slope_fdom*tmp$cdom - C
  }
  
  data.frame(depth = tmp$depth, fluo_cor = fluo_cor, 
             slope_fdom = rep(slope_fdom, length(fluo_cor)),
             C = rep(C, length(fluo_cor)), juld = tmp$juld[1], platform = tmp$platform[1], year = tmp$year[1])
})

### ADD INFO about FDOM-correction for the floats equipped with CDOM ### ===> keep this when writing the article
#fdom_info <- FDOM_OR_MINIMUM_OFFSET_CORRECTION[FDOM_OR_MINIMUM_OFFSET_CORRECTION$platform %in% c(6900807, 6901866),]
#fdom_info <- transform(fdom_info,id=as.numeric(factor(juld)))
# fdom_info <- ddply(fdom_info, ~juld, summarize,
#                    C = mean(C),
#                    slope = mean(slope_fdom))
# 
# meanSlope <- mean(fdom_info$slope)
# sdSlope <- sd(fdom_info$slope)

### END info FDOM

#REPLACE CORRECTED VALUES
profiles[,3] <- FDOM_OR_MINIMUM_OFFSET_CORRECTION$fluo_cor

#QUENCHING CORRECTION
#NOTE : This function was given to me by Marin Cornec (PhD Student at the LOV)
quenching_correction <- function(fluo,depth,MLD) {
  if(is.na(MLD) == FALSE){
    f <- fluo[!is.na(fluo) & depth <= MLD]
    d <- depth[!is.na(fluo) & depth <= MLD]
    zMax <- d[which.max(f)]
    Max <- max(f)
    Corfluo <- fluo
    #Criteria from Schmechtig et al. 2014
    if(!is.na(MLD) & min(f[d<=zMax])<=(0.9*Max)) Corfluo[depth<=zMax] <- Max
    return(Corfluo)
  }else{
    return(fluo)
  }
}

#NPQ correction adapted to daily and nightly (DN) profiles
NPQ_correction_DN_adapted <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id == i,]
  
  if(tmp$day_night[1] == "day"){
    correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[i]) # NPQ correction only applies for daily profiles
  }else{correction <- tmp$fluo}#for nightly profiles
  
  data.frame(fluo_NPQ = correction)
})

profiles$fluo <- NPQ_correction_DN_adapted$fluo_NPQ

#REPLACE < 0 of CHLA by 0 (fair assumption.. ?) ===> MUST BE DONE AFTER THE FDOM CORRECTION (in that case we need to RENORMALIZE CHLA DATA IN PROFILES)
profiles$fluo[which(profiles$fluo < 0)] <- 0

# keep density_profiles and profiles with dcm or not (so all profiles after QC) => 'initial' data
init_profiles <- profiles
init_density_profiles <- clean_remove(density_profiles, unique(density_profiles$juld))


# keep MLD data for 2014-2018 before the fit
init_MLD <- clean_remove(MLDdf, unique(MLDdf$juld))

################################################
### DOWNWELLING PAR  EXTRACTION AND CLEANING ###
################################################

PAR <- ldply(as.list(c("6901866_Mprof.nc", "7900591_Mprof.nc", "7900592_Mprof.nc")),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = F, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  
  FloatInfo <- list(N_PROF=N_PROF,
                    N_LEVELS=N_LEVELS,
                    juld=juld,
                    pres=pres,
                    lon=lon,
                    lat=lat)
  
  pardf <- ExtractVar("DOWNWELLING_PAR", FloatInfo)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  data.frame(depth    = pardf$depth,
             juld     = pardf$juld,
             par     = pardf$value,
             qc       = pardf$qc,
             day      = month.day.year(pardf$juld,c(1,1,1950))$day,
             month    = month.day.year(pardf$juld,c(1,1,1950))$month,
             year     = month.day.year(pardf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(pardf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = pardf$lon,
             lat      = pardf$lat,
             dir      = pardf$dir,
             platform = as.numeric(unique(id)),
             type     = "Argo")
})

PAR <- PAR[-(which(PAR$qc == 4 | PAR$qc == 3)),]

initPAR <- clean_remove(PAR, unique(PAR$juld))

# Due to a visual QC for PAR profiles
initPAR <- initPAR[-which(initPAR$id %in% c(1,3,5,7,8,93,161,175,176,214,318,320,400,495,503,525,527,536,548)),]
initPAR <- clean_remove(initPAR, unique(initPAR$juld))
#REMARK: Pas certain que PAR_0 soit toujours > PAR_NOT_0 (below the surface)

PAR <- clean_remove(initPAR, unique(profiles$juld))

####################################
### DOXY EXTRACTION AND CLEANING ###
####################################

doxy <- ldply(as.list(c("6900807_Mprof.nc", "6901866_Mprof.nc", "7900591_Mprof.nc", "7900592_Mprof.nc")),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = F, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  
  FloatInfo <- list(N_PROF=N_PROF,
                    N_LEVELS=N_LEVELS,
                    juld=juld,
                    pres=pres,
                    lon=lon,
                    lat=lat)
  
  doxydf <- ExtractVar("DOXY", FloatInfo)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  data.frame(depth    = doxydf$depth,
             juld     = doxydf$juld,
             doxy     = doxydf$value,
             qc       = doxydf$qc,
             day      = month.day.year(doxydf$juld,c(1,1,1950))$day,
             month    = month.day.year(doxydf$juld,c(1,1,1950))$month,
             year     = month.day.year(doxydf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(doxydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = doxydf$lon,
             lat      = doxydf$lat,
             dir      = doxydf$dir,
             platform = as.numeric(unique(id)))
})

doxy <- doxy[-(which(doxy$qc == 4 | doxy$qc == 3)),]

initDOXY <- doxy
doxy <- clean_remove(doxy, unique(profiles$juld))

#clean the workspace from some unused variables
rm(argodf, FDOM_OR_MINIMUM_OFFSET_CORRECTION, ncfile, NPQ_correction_DN_adapted, bad_mld_juld, depth_ref, sigma_criteria, filename)
