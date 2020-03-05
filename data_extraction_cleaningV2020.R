#
# This function extracts the data and do some quality controls (i.e. correct chlorophyll a (CHLA) for Non-Photochemical Quenching and for the deep red fluorescence)
# It also derives other parameters such as dissolved oxygen (DOXY) and Photosynthetic Available Radiation (PAR)
#
# In outputs you will find:
# density_profiles* : the density profiles related to the depth of chla profiles (density (sigma) interpolated on depth levels)
# profiles* : the CHLA - CDOM - BBP profiles with some metadata related to it
# MLDdf* : the computed mixed layer depth
# init_density*, init_MLD* and init_profiles* : to save the data (init-) before going further in the script (see the fitting script)
# The * means that they all have the same number of profiles (= 977 after Quality Control)
#
# initDOXY and initPAR are the quality controlled DOXY and PAR profiles
# doxy and PAR are the above dataframes that have been "cross"-checked with the CHLA profiles (based on julian days). Meaning that if a PAR profile has a julian day that is 
# not in the profiles dataframe after quality control, it is removed.

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
library(nortest)

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
DayOrNight <- function(dataframe){
  
  d <- ldply(as.list(unique(dataframe$juld)), function(juld){
    
    print(juld)
    
    origin <- ymd_hms("1950-01-01 00:00:00") # julian day origin for the Argo program
    
    # get temporary profile and metadata
    tmp <- dataframe[dataframe$juld == juld,]
    lat <- tmp$lat[1]
    lon <- tmp$lon[1]
    n <- length(tmp$depth) #number of points in the profile
    
    # determine day or night
    tmp_juld <- origin + juld * 3600 * 24 #conversion of juld in working day and hour
    date <- format(tmp_juld, tz="Europe/Sofia", usetz=TRUE)
    date2 <- as.Date(as.POSIXct(tmp_juld, 'Europe/Sofia')) # convert juld into a a "Date" object
    
    hour <- format(as.POSIXct(strftime(tmp_juld,"%Y-%m-%d %H:%M:%S",tz="Europe/Sofia", usetz = T)),format = "%H")
    data <- data.frame(date = date2, lat = lat, lon = lon) # create a dataframe for the function getSunlightTimes
    sun_cycle <- getSunlightTimes(data = data, keep = c("sunset", "sunrise"), tz = "Europe/Sofia")# get the sunset and sunrise hours
    
    # local sunrise/sunset hours
    hour_sunrise <- as.numeric(format(sun_cycle$sunrise, "%H"))
    hour_sunset <- as.numeric(format(sun_cycle$sunset, "%H"))
    
    if(is.na(hour_sunrise)){ # needed for profiles without latitude/longitude
      sun <- "NA"
    }else if(as.numeric(hour) < hour_sunset & as.numeric(hour) > hour_sunrise){
      sun <- "day"
    }else{
      sun <- "night"
    }
    
    date <- rep(date, times = n)
    sun <- rep(sun, times = n)
    
    data.frame(date = date, day_night = sun)
    
  })
  
  return(d)
}

#4 Function to reinitialize dataframes after cleaning (row numbers, new IDs)
# clean_remove <- function(dataframe, vector_filter){
#   dataframe <- dataframe[dataframe$juld %in% vector_filter,]
#   rownames(dataframe) <- NULL
#   dataframe <- transform(dataframe, id=as.numeric(factor(juld)))
#   dataframe <- dataframe[order(dataframe$id),]
# }

#5 Function to normalize data 
normalize <- function(data){
  data <- (data - min(data, na.rm = T))/(max(data, na.rm = T) - min(data, na.rm = T))
}

#6 Function to QC PAR data (Louis Terrats - Villefranche Oceanographic Laboratory - related article https://journals.ametsoc.org/doi/full/10.1175/JTECH-D-15-0193.1)
radiometric_data_correction <- function(df, parameter) {
  
  ### Above all, define the r² threshold depending on the parameter 
  ### In first, we gather all thresholds defined by Organelli E. et al. 2016
  Thresholds <- setNames(data.frame(matrix(ncol = 3, nrow = 4)), c("Variable", "Threshold 1", "Threshold 2"))
  Thresholds[1,] <- c("ED380", 0.997, 0.999)
  Thresholds[2,] <- c("ED412", 0.997, 0.998)
  Thresholds[3,] <- c("ED490", 0.996, 0.998)
  Thresholds[4,] <- c("PAR", 0.996, 0.998)
  
  ### Then, we keep only thresholds related to the parameter of interest
  Thresholds <- Thresholds[which(Thresholds$Variable == parameter),]
  
  ### If PARAM data are available
  if (all(is.na(df[,parameter])) == FALSE) {
    
    ### Select relevant columns for the PARAM quality checking
    data <- dplyr::select(df, parameter, DEPTH, TIME, LAT, LON) 
    data <- data[which(!is.na(data[,parameter])),]
    
    ### Give default name to the parameter column for writing the function with ease
    colnames(data)[1] <- "PARAM"
    
    ### Compute solar elevation
    solar_elevation <-  sunAngle(as.POSIXct(strptime(unique(data$TIME), "%Y-%m-%d %H:%M:%S"), tz = "UTC"), 
                                 longitude = unique(data$LON), latitude = unique(data$LAT), useRefraction = TRUE)
    
    ### Assign a default flag (= 1)
    data$FLAG <- 1
    data$FEEDBACK <- NA
    
    ### If solar elevation is less than 2 
    if (abs(solar_elevation$altitude) < 2) {
      
      ### All profile flags are equal to 3
      data$FLAG <- 3 
      data$FEEDBACK <- "sun elevation < 2°"
    }
    
    # ### Removing negative spikes (This step is not include in the Emanuele's process)
    # median_values <- decmedian(data$PARAM, order = 2, times = 1, ends = "fill") %>% pastecs::extract(component = "filtered")
    # residuals <- data$PARAM - median_values
    # threshold <- quantile(residuals, probs = 0.10) * 2
    # data$FLAG <- ifelse(residuals < threshold,3,1) 
    # data$FEEDBACK <- ifelse(residuals < threshold,"negative spikes","") 
    
    ### If all profile flags are not equal to 3
    ### Dark value identification
    if ({all(data$FLAG == 3) == FALSE}) {
      
      ### Get the lastest processed data frame 
      dark_values <- data[which(data$FLAG < 3),]
      
      ### If there are at least 5 good quality data
      if ({nrow(dark_values) >= 5}) {
        
        ### Compute the normalized distribution test (Lilliefors test)
        test_results <- lillie.test(dark_values$PARAM) 
        p.value <- test_results$p.value
        
        ### As long as data have not a significant normalized distribution (and sample size is greater than 4, which is one of lillie.test conditions)
        while ({p.value < 0.01} && {nrow(dark_values) > 5}) {
          
          ### Remove the shallowest measure
          dark_values <- dark_values[-1,]
          
          ### And recompure the normalized distribution test
          if (all(dark_values$PARAM == 0)) {
            p.value = 1 
          } else {
            test_results <- lillie.test(dark_values$PARAM) 
            p.value <- test_results$p.value }  
          
        }
        
        ### If dark_values row number is 5, clear out dark_values. 
        ### We will accept the last 4 rows as no dark values (while loop conditions is nrow(dark values) > 5 because of lillie.test requirement which is at least 5 values.
        if (nrow(dark_values) == 5) {
          dark_values <- dark_values[0,]
        }
        
        
        ### Assign 3 to dark value flags
        data$FLAG[which(data$DEPTH %in% dark_values$DEPTH)] <- 3 
        data$FEEDBACK[which(data$DEPTH %in% dark_values$DEPTH)] <- "dark value" 
        
      } else {data$FLAG <- 3} # In case nrow data is < 5
    }
    
    ### If non dark values are less than 5, the profile flag will be 3
    if (nrow(data[which(data$FLAG < 3),]) < 5) {
      data$FLAG <- 3
      data$FEEDBACK <- "not enough data after dark test"
    }
    
    ### If all profile flags are not equal to 3
    if (all(data$FLAG == 3) == FALSE) {
      
      ### Compute the Naperian values of PARAM 
      data$PARAM_ln <- NA
      data$PARAM_ln[which(data$FLAG < 3)] <- suppressWarnings(log(data$PARAM[which(data$FLAG < 3)]))
      
      ### Assign 3 to Na and Inf values created by Napierian logarithm transformation
      data$FEEDBACK[which(data$FLAG < 3 & {is.na(data$PARAM_ln) | is.infinite(data$PARAM_ln)})] <-  "Na or Inf value after log transformation"
      data$FLAG[which(data$FLAG < 3 & {is.na(data$PARAM_ln) | is.infinite(data$PARAM_ln)})] <- 3 
      
      if (nrow(data[which(data$FLAG < 3),]) > 0) {
        
        ### Compute 4 order polynomial regression with non 3 flag values
        a <- lm(data = data[which(data$FLAG < 3),], data$PARAM_ln[which(data$FLAG < 3)] ~ poly(data$DEPTH[which(data$FLAG < 3)], 4, raw = TRUE)) %>% summary() 
        
        ### If r² is less than 0.995
        if (a$r.squared < 0.995) {
          
          ### The profile flag is 3
          data$FLAG <- 3
          data$FEEDBACK <- "r² < 0.995"
          
        } else { # If the r² is more than 0.995
          
          ### Compute the standard deviation as well as mean of residuals from the regression
          residuals_sd <- sd(a$residuals)
          residuals_mean <- mean(a$residuals)
          
          ### If the residual value is out of [mean - 2 sd ; mean + 2 sd], assign 3 to the value's flag
          data$FEEDBACK[which(data$FLAG < 3)] <-  ifelse(abs(a$residuals) > (residuals_mean + 2*residuals_sd), "Residuals out of range (after r² < 0.995 filtering)", "")
          data$FLAG[which(data$FLAG < 3)] <-  ifelse(abs(a$residuals) > (residuals_mean + 2*residuals_sd), 3, 1)
          
          if (nrow(data[which(data$FLAG < 3),]) > 0) {
            
            ### Then, compute a second 4 order polynomial regression with non 3 flag values
            b <- lm(data = data[which(data$FLAG < 3),], data$PARAM_ln[which(data$FLAG < 3)] ~ poly(data$DEPTH[which(data$FLAG < 3)], 4, raw=TRUE)) %>% summary() 
            
            ### If r² is less than the threshold 1
            if (round(b$r.squared, digits = 3) < as.numeric(Thresholds[,2])) {
              
              ### The profile flag will be 3
              data$FLAG <- 3
              data$FEEDBACK <- "r² < Threshold 1"
              
            } else { ### or if the r² is more than the threshold 1
              
              ### If the r² is less than the threshold 2 (= in the range [Threshold 1 ; Threshold 2])
              if (round(b$r.squared, digits = 3) < as.numeric(Thresholds[,3])) {
                
                ### Compute the mean as well as standard deviation of residuals from the second regression
                residuals_sd <- sd(b$residuals)
                residuals_mean <- mean(b$residuals)
                
                ### If the residual value is out of the range  [mean - 2 sd ; mean + 2 sd], the value's flag will be 3 
                ### If not, the value's flag will be 2
                data$FEEDBACK[which(data$FLAG < 3)] <-  ifelse(abs(b$residuals) > (residuals_mean + 2*residuals_sd), "Residuals out of range (after r² > threshold 1 filtering)", "")
                data$FLAG[which(data$FLAG < 3)] <-  ifelse(abs(b$residuals) > (residuals_mean + 2*residuals_sd), 3, 2)
                
              } 
              
              ### If the r² is more than the threshold 2
              if (round(b$r.squared, digits = 3) >= as.numeric(Thresholds[,3])) { 
                
                ### If the residual's value is less than [mean - 1 sd ; mean + 1 sd], keep 1 to the value's flag
                ### If the residual's value is more than [mean + 2 sd; mean - 2 sd], give 3 to the value's flag 
                ### Otherwise assign 2 to the value's flag
                data$FEEDBACK[which(data$FLAG < 3)] <-  ifelse(abs(b$residuals) > residuals_mean + 2*residuals_sd, "Residuals out of range (after r² > threshold 2 filtering)", 
                                                               ifelse(abs(b$residuals) > residuals_mean + 1*residuals_sd, "", ""))
                data$FLAG[which(data$FLAG < 3)] <-  ifelse(abs(b$residuals) > residuals_mean + 2*residuals_sd, 3, 
                                                           ifelse(abs(b$residuals) > residuals_mean + 1*residuals_sd, 2, 1))
                
              }
              
            }
          }
        }} 
      
      colnames(data)[which(names(data) == "PARAM")] <- parameter
      colnames(data)[which(names(data) == "FLAG")] <- paste(parameter, "_FLAG", sep = "")
      df <- full_join(df, data[,c("DEPTH", parameter, paste(parameter, "_FLAG", sep = ""), "FEEDBACK")], by = c("DEPTH", parameter))
      
      return(df)
    } }
}

##########
# FLOATS #
##########

filename <- c("6900807_Mprof2020.nc","6901866_Mprof2020.nc","7900591_Mprof2020.nc","7900592_Mprof2020.nc","6903240_Mprof2020.nc")

#####################
# GET PROFILES DATA #
#####################

# Acquisition of CHLA and BBP (non adjusted !) and CDOM (if available) profiles from each float
profiles <- ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = F, suppress_dimvals = FALSE)
  
  print(file)
  
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
  
  if (file == "7900591_Mprof2020.nc" | file == "7900592_Mprof2020.nc"){ #change M for merged or S for synthetic profiles
    cdomdf <- as.data.frame(rep(NA, length(chladf$value)))
    colnames(cdomdf) <- "value"
  }else{
    cdomdf <- ExtractVar("CDOM", FloatInfo)
  }
  
  # for a dimension mismatch issue occurring with FLOAT WMO 6900807 (between BBP and CHLA)
  if(file == "6900807_Mprof2020.nc"){
    chladf <- join(bbpdf, chladf, by = c("alevel", "aprofile"))
    chladf <- chladf[,11:18]
    cdomdf <- join(bbpdf, cdomdf, by = c("alevel", "aprofile"))
    cdomdf <- cdomdf[,11:18]
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  day_night <- DayOrNight(chladf)
  
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
             day_night = day_night$day_night,
             date_long = day_night$date,
             platform = as.numeric(unique(id))
  )
})

# save all profiles
#save.image(file = 'allprofilesbeforeQC.Rdata')
#load(file = 'allprofilesbeforeQC.Rdata')

######################
### DATA CLEANING  ###
######################

# REMOVE DATA OF 2013 and 2020 (incomplete years)
profiles <- profiles[!(profiles$year == 2013 | profiles$year == 2020),]

# REMOVE DESCENT PROFILES
profiles <- profiles[-(which(profiles$dir == "D")),]

# REMOVE PROFILES WITH NO LAT/LON
profiles <- profiles[-(which(is.na(profiles$lat == T))),]

# Get some basic info on each profile
argodf <- ddply(profiles,~juld,summarize,
                depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],# to not take into account surface depth where chlorophyll can be minimum
                bottomdepth = max(depth),
                min_depth = min(depth),
                platform = platform[1],
                month = month[1],
                year = year[1], DOY = DOY[1],
                n = length(fluo)) #n: number of points in the profile

# Get profiles id that do not answer to some criteria
crit1 <- argodf[(argodf$min_depth > 5),]# all stuck profiles except one (add QC on pressure to get rid of them directly)
crit2 <- argodf[(argodf$bottomdepth < 40),] # 50m works as well, justification: le dcm peut être en-dessous de 50m 

# remove a profile (.6 24895.96 0.8784 0.00916 2.7240 2 2 2 2018 59 34.85083 44.71711 A night 6900807 652) that shows a DCM at > 150m deep (cannot correct with FDOM method) -> remove it
crit3 <- argodf[argodf$platform == 6900807 & argodf$DOY == 59 & argodf$year == 2018,]

criteriadf <- rbind(crit1, crit2, crit3)
criteriadf <- criteriadf[(!duplicated(criteriadf)),] # this line can be removed because it does not change anything for our case

#remove profiles that have been categorized as bad profiles
profiles <- profiles[-which(profiles$juld %in% criteriadf$juld),]

# clean workspace
rm(crit1, crit2, crit3, criteriadf)

# REMOVE BAD DATA (I.E. QC = 4 FOR NEGATIVE SPIKES, JUMPS, ETC.)
profiles <- profiles[-(which(profiles$qc_fluo == 4)),] # bad chla data
profiles <- profiles[-(which(profiles$qc_bbp %in% c(3,4))),] # bad bbp data

# MEDIAN FILTER TO REMOVE POSITIVE SPIKES (THEY MAY RETAIN INFORMATION BUT STRONG SPIKES WILL NOT BE RETAINED FOR THIS STUDY)
smoothed_fluo <- ldply(as.list(unique(profiles$juld)), function(i){
  tmp <- profiles[profiles$juld == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

smoothed_bbp <- ldply(as.list(unique(profiles$juld)), function(i){
  tmp <- profiles[profiles$juld == i,]
  tmp <- mmed(tmp$bbp, 5)
  data.frame(smoothed_bbp = tmp)
})

smoothed_cdom <- ldply(as.list(unique(profiles$juld)), function(i){
  tmp <- profiles[profiles$juld == i,]
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
rm(smoothed_bbp, smoothed_fluo, smoothed_cdom, ncfile, argodf)

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

# clean data
density_profiles <- filter(density_profiles, juld %in% profiles$juld)

#COMPUTE MLD (criteria can thus be modified according to the definition of MLD that you choose)
sigma_criteria <- 0.03
depth_ref <- 10

MLDdf <- ldply(as.list(unique(density_profiles$juld)), function(i){
  
  tmp <- density_profiles[density_profiles$juld == i,]
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
    sigmaMaxMLD <- mean(tmp$sigma[which(tmp$depth == MLD)], na.rm =T) #check this, pq max mld? #####################################################################
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

# clean data
density_profiles <- filter(density_profiles, juld %in% profiles$juld)
MLDdf <- filter(MLDdf, juld %in% profiles$juld)

# add a column to the 'profiles' dataframe with OLD_CHLA (will save the uncorrected chla, can be useful for some plots)
profiles$old_fluo <- profiles$fluo

## In reality, it could be a script based not on the ID but on the JULD (avoid all the rownames init and so forth.. To meditate)

# FDOM-BASED CORRECTION AND MINIMUM-OFFSET CORRECTION
FDOM_OR_MINIMUM_OFFSET_CORRECTION <- ldply(as.list(unique(profiles$juld)), function(i){
  tmp <- profiles[profiles$juld == i,]
  MaxDepth <- max(tmp$depth) # Max depth of the profile
  TopDepth <- tmp$depth[which.max(tmp$fluo):length(tmp$fluo)][which.min(tmp$fluo[which.max(tmp$fluo):length(tmp$fluo)])]# Apparent minimum 
  
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
# fdom_info <- FDOM_OR_MINIMUM_OFFSET_CORRECTION[FDOM_OR_MINIMUM_OFFSET_CORRECTION$platform %in% c(6900807, 6901866),]
# fdom_info <- transform(fdom_info,id=as.numeric(factor(juld)))
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

#NPQ correction adapted to daytime and nighttime (DN) profiles
NPQ_correction_DN_adapted <- ldply(as.list(unique(profiles$juld)), function(i){
  
  tmp <- profiles[profiles$juld == i,]
  
  if(tmp$day_night[1] == "day"){
    correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[which(MLDdf$juld == i)]) # NPQ correction only applies for daily profiles
  }else{correction <- tmp$fluo}#for nighttime profiles
  
  data.frame(fluo_NPQ = correction)
})

profiles$fluo <- NPQ_correction_DN_adapted$fluo_NPQ

#REPLACE < 0 of CHLA by 0 (fair assumption.. ?) ===> MUST BE DONE AFTER THE FDOM CORRECTION (in that case we need to RENORMALIZE CHLA DATA IN PROFILES)
profiles$fluo[which(profiles$fluo < 0)] <- 0

# keep density_profiles and profiles with dcm or not (so all profiles after QC) => 'initial' data
init_profiles <- profiles
init_density_profiles <- filter(density_profiles, juld %in% profiles$juld)

# keep MLD data for 2014-2018 before the fit
init_MLD <- filter(MLDdf, juld %in% profiles$juld)

################################################
### DOWNWELLING PAR  EXTRACTION AND CLEANING ###
################################################

# 6900807 does not have PAR
PAR <- ldply(as.list(c("6901866_Mprof2020.nc","7900591_Mprof2020.nc","7900592_Mprof2020.nc","6903240_Mprof2020.nc")),function(file){

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

# QUALITY CONTROL (THANKS TO LOUIS TERRATS'S SCRIPT, see also the article related to his function)
colnames(PAR) <- c("DEPTH","juld","PAR","QC","DAY","MONTH","YEAR","DOY","LON","LAT","DIR","PLATFORM","TYPE")

origin <- ymd_hms("1950-01-01 00:00:00") # julian day origin for the Argo program
tmp_juld <- origin + PAR$juld * 3600 * 24 #conversion of juld in working day and hour
PAR$TIME <- tmp_juld
#remove lines with lat/lon = NAN 
PAR <- PAR[which(!is.na(PAR$LAT)),]  #for my case, when LAT is not present, so is LON
PAR <- PAR[PAR$YEAR %in% c(2014:2019),]
PAR <- transform(PAR, id=as.numeric(factor(TIME)))

#### TAKES WAY TO MUCH TIME TO COMPUTE WHEN WE GIVE ALL PROFILES AT ONCE (BATCH)

# Split it? Maybe it's too much for the memory?
init_dataframe <- radiometric_data_correction(PAR[PAR$id == 1,], "PAR")

for (i in 2:length(unique(PAR$id))){
  print(i)
  tmp <- radiometric_data_correction(PAR[PAR$id == i,], "PAR")
  init_dataframe <- rbind(init_dataframe, tmp)
}

final_PAR <- init_dataframe
initPAR <- final_PAR[final_PAR$PAR_FLAG %in% c(1,2),]

colnames(initPAR) <- c("depth","juld","par","qc","day","month","year","DOY","lon","lat","dir","platform","type","time","id","par_flag","feedback")

####################################
### DOXY EXTRACTION AND CLEANING ###
####################################

doxy <- ldply(as.list(c("6901866_Mprof2020.nc","7900591_Mprof2020.nc","7900592_Mprof2020.nc","6903240_Mprof2020.nc")),function(file){
  
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

#clean the workspace from some unused variables
rm(FDOM_OR_MINIMUM_OFFSET_CORRECTION, ncfile, NPQ_correction_DN_adapted, 
   bad_mld_juld, depth_ref, sigma_criteria, filename, i, origin, tmp, tmp_juld, final_PAR, init_dataframe)

##########################################
# Now that the data are QCed and cleaned #
##########################################

###################################### This is needed because (T,S) are from C-files and BGC data are from B-files 
# DENSITY PROFILES WITH CHLA AND BBP # (merged in the M-files but this may be solvable with the S-files?)
######################################

# Let's a dataframe based on potential density anomaly, CHLA and BBP
density_CHLA_BPP_profiles <- ldply(as.list(unique(density_profiles$juld)), function(i){
  # rmk here density and init profiles (also profiles and init profiles are still the same)
  tmp <- density_profiles[density_profiles$juld == i,]
  tmp2 <- profiles[profiles$juld == i,]
  chla_bbp_density <- approx(tmp$depth, tmp$sigma, tmp2$depth)$y #interpolation
  tmp2 <- cbind(tmp2, chla_bbp_density)
  
  data.frame(sigma = tmp2$chla_bbp_density)
})

profiles$sigma <- density_CHLA_BPP_profiles$sigma
init_profiles$sigma <- density_CHLA_BPP_profiles$sigma
rm(density_CHLA_BPP_profiles)


# Normalization of CHLA and BBP 
norm_tmp <- ldply(as.list(unique(profiles$juld)), function(i){
  tmp <- profiles[profiles$juld == i,]
  chla_norm <- normalize(tmp$fluo)
  bbp_norm <- normalize(tmp$bbp)
  sigma_norm <- normalize(tmp$sigma)
  data.frame(chla_norm = chla_norm, bbp_norm = bbp_norm, sigma_norm = sigma_norm)
})

profiles$fluo_norm <- norm_tmp$chla_norm
profiles$bbp_norm <- norm_tmp$bbp_norm
profiles$sigma_norm <- norm_tmp$sigma_norm

init_profiles$fluo_norm <- norm_tmp$chla_norm
init_profiles$bbp_norm <- norm_tmp$bbp_norm
init_profiles$sigma_norm <- norm_tmp$sigma_norm

norm_doxy <- ldply(as.list(unique(initDOXY$juld)), function(i){
  tmp <- initDOXY[initDOXY$juld == i,]
  doxy_norm <- normalize(tmp$doxy)
  data.frame(depth = tmp$depth, doxy_norm = doxy_norm)
})

initDOXY$doxy_norm <- norm_doxy$doxy_norm
doxy <- filter(initDOXY, juld %in% unique(profiles$juld))

norm_PAR <- ldply(as.list(unique(initPAR$juld)), function(i){
  tmp <- initPAR[initPAR$juld == i,]
  par_norm <- normalize(tmp$par)
  data.frame(depth = tmp$depth, par_norm = par_norm)
})

initPAR$par_norm <- norm_PAR$par_norm
PAR <- filter(initPAR, juld %in% unique(profiles$juld))


rm(norm_tmp, norm_doxy, norm_PAR)