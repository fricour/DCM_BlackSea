# Script for extracting data from the new float 6903240

#############
# LIBRARIES #
#############

library(ncdf4)
library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(ggforce)

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

########
# DATA #
########

filename <- "6903240_Mprof.nc"

data <- ldply(as.list(filename),function(file){
  
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

#REMOVE BAD DATA
data <- data[-(which(data$qc_fluo == 4)),] 
rownames(data) <- NULL
data <- data[-(which(data$depth < 0)),]

#Day of deployment
data <- data[data$DOY == 88 & data$year == 2018,]

smoothed_fluo <- ldply(as.list(unique(data$juld)), function(i){
  tmp <- data[data$juld == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

smoothed_cdom <- ldply(as.list(unique(data$juld)), function(i){
  tmp <- data[data$juld == i,]
  tmp <- mmed(tmp$cdom, 5)
  data.frame(smoothed_cdom = tmp)
})

#REPLACE FLUO WITH SMOOTHED FLUO & CDOM
data[,3] <- smoothed_fluo
data[,4] <- smoothed_cdom

data <- clean_remove(data, unique(data$juld))

data_info <- ddply(data,~juld,summarize,
                    qc_fluo = qc_fluo[which.max(fluo)],
                    depthmax = depth[which.max(fluo)],
                    maxvalue = fluo[which.max(fluo)],
                    depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],
                    integration = sum(fluo),
                    bottomdepth = max(depth),
                    min_depth = min(depth),
                    dir = dir[1],
                    lon=mean(lon),
                    lat=mean(lat),
                    platform = platform[1],
                    day = day[1], month = month[1],
                    year = year[1], DOY = DOY[1])

rm(ncfile, smoothed_cdom, smoothed_fluo)

# compute density profiles
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

density_profiles <- clean_remove(density_profiles, unique(data$juld))
data_info <- transform(data_info,id=as.numeric(factor(juld)))
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))

#MLD computation
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

# manual adjust because density inversion
MLDdf$MLD[3] <- 23.5

# add a column to the 'profiles' dataframe with OLD_CHLA (will save the uncorrected chla, can be useful for some plots)
data$old_fluo <- data$fluo

#EACH PROFILES HAS BEEN VALIDATED SO FAR ==> CALIBRATION EXERCISE STARTS NOW
#LET'S DO IT ONLY FOR THE DEEP PROFILE
#FDOM-based method (Xing et al. 2017)

#keep only the validation profile (i.e. id = 1)
data <- data[data$id == 3,]
data_info <- clean_remove(data_info, unique(data$juld))
MLDdf <- clean_remove(MLDdf, unique(data$juld))
data <- transform(data,id=as.numeric(factor(juld)))

#FDOM-BASED CORRECTION AND MINIMUM-OFFSET CORRECTION
FDOM_OR_MINIMUM_OFFSET_CORRECTION <- ldply(as.list(1:length(unique(data$id))), function(i){
  tmp <- data[data$id == i,]
  MaxDepth <- data_info$bottomdepth[i]#Max depth of the profile
  TopDepth <- data_info$depthmin[i]#Apparent minimum 
  
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

#USE THE FDOM
data[,3] <- FDOM_OR_MINIMUM_OFFSET_CORRECTION$fluo_cor

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
NPQ_correction_DN_adapted <- ldply(as.list(1:length(unique(data$id))), function(i){
  
  tmp <- data[data$id == i,]
  
  if(tmp$day_night[1] == "day"){
    correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[i]) # NPQ correction only applies for daily profiles
  }else{correction <- tmp$fluo}#for nightly profiles
  
  data.frame(fluo_NPQ = correction)
})

data$fluo_npq <- NPQ_correction_DN_adapted$fluo_NPQ

#REPLACE < 0 of CHLA by 0 (fair assumption.. ?) ===> MUST BE DONE AFTER THE FDOM CORRECTION (in that case we need to RENORMALIZE CHLA DATA IN PROFILES)
data$fluo[which(data$fluo < 0)] <- 0
data$fluo_npq[which(data$fluo_npq < 0)] <- 0

########################
# Comparison with HPLC #
########################

#HPLC DATA FROM DEPLOYMENT
depth_hplc <- c(0,30,50,70,90,100,140,200,250,500,750,1000)
chla_hplc <- c(0.6111,0.6433,0.3090,0.0675,0.0229,0.0168,0.012,0.0032,
               0.0035, 0.0023, 0.0039, 0.0042)


hplc <- as.data.frame(cbind(depth_hplc, chla_hplc))

deployment <- ggplot(data, aes(x = fluo_npq, y = depth, colour = "FDOM + NPQ")) + geom_path(size = 0.5) +
  scale_y_reverse() + geom_path(data = data, aes(x = fluo, y = depth, colour = "FDOM")) +
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)") + geom_hline(yintercept = data_info$depthmin[1], color = "black", lty = 2) +
  annotate("text", x = .8, y = data_info$depthmin[1]-20, label = "CHLA MIN.") +
  geom_path(data = data, aes(x = old_fluo, y = depth, colour = "RAW")) + 
  geom_point(data=hplc, aes(y = depth_hplc, x = chla_hplc), colour = "red", size = 2, shape = 15)  + 
  geom_hline(yintercept = MLDdf$MLD[1], color = "black", lty = 2) + scale_color_discrete(name = "Method")+
  annotate("text", x = .2, y = MLDdf$MLD[1]-20, label = "MLD") + 
  theme(legend.position = "none")

zoom_deployment <- ggplot(data, aes(x = fluo_npq, y = depth, colour = "FDOM + NPQ")) + geom_path(size = 0.5) +
  scale_y_reverse(limits = c(100,0)) + geom_path(data = data, aes(x = fluo, y = depth, colour = "FDOM")) +
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)") + geom_hline(yintercept = data_info$depthmin[1], color = "black", lty = 2) +
  annotate("text", x = .8, y = data_info$depthmin[1]-2, label = "CHLA MIN.") +
  geom_path(data = data, aes(x = old_fluo, y = depth, colour = "RAW")) + 
  geom_point(data=hplc, aes(y = depth_hplc, x = chla_hplc), colour = "red", size = 2, shape = 15)  + 
  geom_hline(yintercept = MLDdf$MLD[1], color = "black", lty = 2) + scale_color_discrete(name = "Method") +
  annotate("text", x = .2, y = MLDdf$MLD[1]-2, label = "MLD") +
  theme(axis.title.y=element_blank()) 

grid.arrange(deployment, zoom_deployment, nrow = 1, ncol = 2)


# COMPUTE STATS BELOW

indexdepth <- vector()

for (i in 1:length(hplc$chla)){
  indexdepth[i] <- which.min(abs(data$depth - hplc$depth_hplc[i]))
}

# RMSE DEEP LAYER FDOM CORRECTED PROFILE
index <- which(hplc$depth_hplc > data_info$depthmin & hplc$depth_hplc < data_info$bottomdepth)
deep_layer <- indexdepth[index]
rmse_deep_layer <- sqrt(sum((hplc$chla_hplc[index] - data$fluo[deep_layer])^2)/length(hplc$depth_hplc[index]))

# SAME FOR THE RAW PROFILE
rmse_deep_layer_old <- sqrt(sum((hplc$chla_hplc[index] - data$old_fluo[deep_layer])^2)/length(hplc$depth_hplc[index]))

############# 
# OLD STUFF #
#############

# rmse1 <- data$fluo[indexdepth]
# rmseFDOM1 <- validationdf$fluo_FDOM[indexdepth]
# rmseFDOMNPQ1 <- validationdf$fluo_FDOM_NPQ[indexdepth]
# 
# #########
# 
# #RMSE
# rmse <- sqrt(sum((hplc$chla_hplc - rmse1)^2)/length(hplc$depth_hplc))
# #RMSE for top layer correction
# rmseTOP <- sqrt(sum((hplc$chla_hplc[1:5] - rmse1[1:5])^2)/length(hplc$depth_hplc[1:5]))
# #for deep
# rmseDEEP <- sqrt(sum((hplc$chla_hplc[6:length(depth_hplc)] - rmse1[6:length(depth_hplc)])^2)/length(hplc$depth_hplc[6:length(depth_hplc)]))
# #below 50m
# #rmse <- sqrt(sum((hplc$chla_hplc[4:9] - rmse_chla_vec[4:9])^2)/length(hplc$depth_hplc[4:9]))
# 
# #RMSE with FDOM
# rmseFDOM <- sqrt(sum((hplc$chla_hplc - rmseFDOM1)^2)/length(hplc$depth_hplc))
# #RMSE for top layer correction
# rmseFDOMTOP <- sqrt(sum((hplc$chla_hplc[1:5] - rmseFDOM1[1:5])^2)/length(hplc$depth_hplc[1:5]))
# #for deep
# rmseFDOMDEEP <- sqrt(sum((hplc$chla_hplc[6:length(depth_hplc)] - rmseFDOM1[6:length(depth_hplc)])^2)/length(hplc$depth_hplc[6:length(depth_hplc)]))
# 
# #RMSE with FDOM AND NPQ (note que ça ne changera que la valeur TOP and ALL (not DEEP))
# rmseFDOMNPQ <- sqrt(sum((hplc$chla_hplc - rmseFDOMNPQ1)^2)/length(hplc$depth_hplc))
# #RMSE for top layer correction
# rmseFDOMNPQTOP <- sqrt(sum((hplc$chla_hplc[1:5] - rmseFDOMNPQ1[1:5])^2)/length(hplc$depth_hplc[1:5]))
# #for deep
# rmseFDOMNPQDEEP <- sqrt(sum((hplc$chla_hplc[6:length(depth_hplc)] - rmseFDOMNPQ1[6:length(depth_hplc)])^2)/length(hplc$depth_hplc[6:length(depth_hplc)]))
# 
# 
# # Xing2017 - Fig 4a
# C <- FDOMdf_validation$C[choice]
# slope <- FDOMdf_validation$slope_fdom[choice]
# min_apparent <- FDOMdf_validation$min_apparent[choice]
# data_below_chla_min <- validation_profile[validation_profile$depth >= min_apparent,]
# ggplot(validation_profile, aes(x = cdom, y = fluo)) + geom_point() + geom_abline(slope = slope, intercept = C) + geom_point(data = data_below_chla_min, colour = "red")
# 
# library(ggpubr)
# ggscatter(data_below_chla_min, x = "cdom", y = "fluo", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")
