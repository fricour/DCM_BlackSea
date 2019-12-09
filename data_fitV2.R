# Mathematical function for CHLA profiles
# Following Mignot et al. 2011 and Carranza et al. 2018
fgauss_expo <- function(z, Fsurf, Zdemi, Fmax, Zmax, dz){
  Fsurf*exp((-log(2)/Zdemi)*z) + Fmax*exp(-(z-Zmax)^2/dz^2)
}

fsigmoid <- function(z, Fsurf, Zdemi, s){
  Fsurf*(1/(1+exp((Zdemi-z)*s)))
}

fexpo <- function(z, Fsurf, Zdemi){
  Fsurf*exp((-log(2)/Zdemi)*z)
}

fgauss <- function(z, Fmax, Zmax, dz){
  Fmax*exp(-(z-Zmax)^2/dz^2)
}

fgauss_sigmoid <- function(z, Fsurf, Zdemi, s, Fmax, Zmax, dz){
  Fmax*exp(-(z-Zmax)^2/dz^2) + Fsurf*(1/(1+exp((Zdemi-z)*s)))
}

#########################
##### BEGIN THE FIT #####
#########################

fit <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id==i,]
  
  #data needed for R^2 and R^2 adjusted
  mean_fluo <- mean(tmp$fluo, na.rm=T) #mean fluo over all profile
  ss_tot <- sum((tmp$fluo - mean_fluo)^2, na.rm=T) #for R^2 computation
  n <- length(tmp$depth) #number of data points
  #Arthur's suggestion => Fsurf should be a mean of the fluo in the MLD
  index_MLD <- which.min(tmp$depth <= MLDdf$MLD[i]) 
  
  ##############
  # PARAMETERS #
  ##############
  
  #Parameters estimation
  Fsurf <- mean(tmp$fluo[1:index_MLD], na.rm = T)
  Zdemi <- tmp$depth[which.max(tmp$fluo <=  Fsurf/2)]
  if(MLDdf$MLD[i] > Zdemi){
    Zdemi <- MLDdf$MLD[i]
  }
  maxindex <- which.max(tmp$fluo)
  Zmax <- tmp$depth[maxindex]
  Fmax <- tmp$fluo[maxindex]
  dz <- 5 # First guess, cannot do better
  s <- -0.01 #guess without explication like dz, it is fixed and adapted if initial values do not work
 
  #############################
  # GAUSSIAN-EXPONENTIAL (GE) #
  #############################
  
  res_GE <- tryCatch(nlsLM(fluo ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  
  if(inherits(res_GE, "error")){# handle error if initial dz is not good
    dz <- 10
    res_GE <- tryCatch(nlsLM(fluo ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  #R^2 and R^2 adjusted
  ge <- fgauss_expo(tmp$depth,coef(res_GE)["Fsurf"], coef(res_GE)["Zdemi"], coef(res_GE)["Fmax"], coef(res_GE)["Zmax"], coef(res_GE)["dz"])
  ss_res_ge <- sum((ge - tmp$fluo)^2, na.rm=T)  
  rcoef_ge <- 1-(ss_res_ge/ss_tot)
  rcoef_ge <- 1 - (1-rcoef_ge)*(n-1)/(n-5-1) # final R^2 ADJUSTED
  
  ###############
  # SIGMOID (S) #
  ###############
  
  res_S <- tryCatch(nlsLM(fluo ~ fsigmoid(depth, Fsurf, Zdemi, s),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  if(inherits(res_S, "error")){#added due to a bug
    Zdemi <- 30
    #s <- -0.01
    res_S <- tryCatch(nlsLM(fluo ~ fsigmoid(depth, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  
  sigmoid <- fsigmoid(tmp$depth, coef(res_S)["Fsurf"], coef(res_S)["Zdemi"], coef(res_S)["s"])
  ss_res_sigmoid <- sum((sigmoid - tmp$fluo)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)               
  rcoef_sigmoid <- 1 - (1-rcoef_sigmoid)*(n-1)/(n-3-1)
  
  ################
  # GAUSSIAN (G) #
  ################
  
  res_G <- tryCatch(nlsLM(fluo ~ fgauss(depth, Fmax, Zmax, dz),
                        start = c(Fmax = Fmax, Zmax = Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  gauss <- fgauss(tmp$depth, coef(res_G)["Fmax"], coef(res_G)["Zmax"], coef(res_G)["dz"])
  ss_res_gauss <- sum((gauss - tmp$fluo)^2, na.rm=T)                
  rcoef_gauss <- 1-(ss_res_gauss/ss_tot)  
  rcoef_gauss <- 1 - (1-rcoef_gauss)*(n-1)/(n-3-1)
  
  ###################
  # EXPONENTIAL (E) #
  ###################
  
  res_E <- tryCatch(nlsLM(fluo ~ fexpo(depth, Fsurf, Zdemi),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  expo <- fexpo(tmp$depth, coef(res_E)["Fsurf"], coef(res_E)["Zdemi"])
  ss_res_expo <- sum((expo - tmp$fluo)^2, na.rm=T)
  rcoef_expo <- 1-(ss_res_expo/ss_tot)
  rcoef_expo <- 1 - (1-rcoef_expo)*(n-1)/(n-2-1)
  
  #########################
  # GAUSSIAN-SIGMOID (GS) #
  #########################
  
  res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                           Fmax, Zmax, dz),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  if(inherits(res_GS, "error")){
    Zdemi <- 30#added due to a bug
    s <- -0.1 #same as before, guess without explainations
    res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                             Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res_GS, "error")){
    dz <- 10
    s <- -0.1
    res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                             Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  gs <- fgauss_sigmoid(tmp$depth, coef(res_GS)["Fsurf"], coef(res_GS)["Zdemi"], coef(res_GS)["s"], coef(res_GS)["Fmax"], coef(res_GS)["Zmax"], coef(res_GS)["dz"])
  ss_res_gs <- sum((gs - tmp$fluo)^2, na.rm=T)
  rcoef_gs <- 1-(ss_res_gs/ss_tot)
  rcoef_gs <- 1 - (1-rcoef_gs)*(n-1)/(n-6-1)
  
  #######################
  # SHAPE DETERMINATION #
  #######################
  
  #which is the best fit (the order was chosen so that if multiple values are equal, the first is taken) ## can be improved?
  bestfit <- which.max(c(rcoef_sigmoid, rcoef_expo, rcoef_gauss, rcoef_ge, rcoef_gs))

  #rejection criteria 
  Rcrit <- 0.9
  
  if(c(rcoef_sigmoid, rcoef_expo, rcoef_gauss, rcoef_ge, rcoef_gs)[bestfit] < Rcrit){
    Zmax <- NA #Zmax is possible but the fit is too bad
    r2adj <- NA
    shape <- "O" # "O" for "Other"
    chla_at_dcm <- NA
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 1){
    Zmax <- NA #no Zmax
    r2adj <- rcoef_sigmoid
    shape <- "S"
    chla_at_dcm <- NA
    #no need for AIC (see Arthur for AIC but we did an R2 *adjusted* instead)
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 2){
    Zmax <- NA
    r2adj <- rcoef_expo
    shape <- "E"
    chla_at_dcm <- NA
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 3){
    Zmax <- coef(res_G)["Zmax"]
    r2adj <- rcoef_gauss
    shape <- "G"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))] 
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- coef(res_G)["Fmax"]
    Zmax <- coef(res_G)["Zmax"]
    dz <- coef(res_G)["dz"]
  }else if(bestfit == 4){
    Zmax <- coef(res_GE)["Zmax"]
    r2adj <- rcoef_ge
    shape <- "GE"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))]
    Fsurf <- coef(res_GE)["Fsurf"]
    Zdemi <- coef(res_GE)["Zdemi"]
    s <- NA
    Fmax <- coef(res_GE)["Fmax"]
    Zmax <- coef(res_GE)["Zmax"]
    dz <- coef(res_GE)["dz"]
  }else{
    Zmax <- coef(res_GS)["Zmax"]
    r2adj <- rcoef_gs
    shape <- "GS"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))]
    Fsurf <- coef(res_GS)["Fsurf"]
    Zdemi <- coef(res_GS)["Zdemi"]
    s <- coef(res_GS)["s"]
    Fmax <- coef(res_GS)["Fmax"]
    Zmax <- coef(res_GS)["Zmax"]
    dz <- coef(res_GS)["dz"]
  }
  
  show(i)
  
  data.frame(shape = shape, r2adj = r2adj, Zmax = Zmax, id = tmp$id[1], juld = tmp$juld[1],
             lon = tmp$lon[1], lat = tmp$lat[1], day = tmp$day[1],
             month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$platform[1], chla_surf = tmp$fluo[1], #chla_surf is the chla value at the surface of the profile
             day_night = tmp$day_night[1], MLD = MLDdf$MLD[i], chla_dcm = chla_at_dcm, Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, dz = dz,
             date = tmp$date_long[1]) 
})

# save fit
init_fit <- fit


###################################### This is needed because (T,S) are from C-files and BGC data are from B-files 
# DENSITY PROFILES WITH CHLA AND BBP # (merged in the M-files but this may be solvable with the S-files?)
######################################

# Let's a dataframe based on potential density anomaly, CHLA and BBP
density_CHLA_BPP_profiles <- ldply(as.list(1:length(unique(density_profiles$id))), function(i){
  
  tmp <- density_profiles[density_profiles$id == i,]
  tmp2 <- profiles[profiles$id == i,]
  chla_bbp_density <- approx(tmp$depth, tmp$sigma, tmp2$depth)$y #interpolation
  tmp2 <- cbind(tmp2, chla_bbp_density)
  
  data.frame(sigma = tmp2$chla_bbp_density)
})

profiles$sigma <- density_CHLA_BPP_profiles$sigma
rm(density_CHLA_BPP_profiles)

# Let's do the same as above for init_profiles (needed for the test of Arthur)
init_density_profiles <- clean_remove(init_density_profiles, unique(init_profiles$juld)) # get the same number of profiles at init ====> I think we can remove this line

density_CHLA_BPP_profiles <- ldply(as.list(1:length(unique(init_density_profiles$id))), function(i){
  
  tmp <- init_density_profiles[init_density_profiles$id == i,]
  tmp2 <- init_profiles[init_profiles$id == i,]
  chla_bbp_density <- approx(tmp$depth, tmp$sigma, tmp2$depth)$y 
  tmp2 <- cbind(tmp2, chla_bbp_density)
  
  data.frame(sigma = tmp2$chla_bbp_density)
})

init_profiles$sigma <- density_CHLA_BPP_profiles$sigma
rm(density_CHLA_BPP_profiles)

# Normalization of CHLA and BBP
norm_tmp <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id == i,]
  chla_norm <- normalize(tmp$fluo)
  bbp_norm <- normalize(tmp$bbp)
  sigma_norm <- normalize(tmp$sigma)
  data.frame(chla_norm = chla_norm, bbp_norm = bbp_norm, sigma_norm = sigma_norm)
})

norm_doxy <- ldply(as.list(1:length(unique(doxy$id))), function(i){
  tmp <- doxy[doxy$id == i,]
  doxy_norm <- normalize(tmp$doxy)
  data.frame(depth = tmp$depth, doxy_norm = doxy_norm)
})

norm_PAR <- ldply(as.list(1:length(unique(PAR$id))), function(i){
  tmp <- PAR[PAR$id == i,]
  par_norm <- normalize(tmp$par)
  data.frame(depth = tmp$depth, par_norm = par_norm)
})

# insert normalized CHLA and BBP into PROFILES
profiles$fluo_norm <- norm_tmp$chla_norm
profiles$bbp_norm <- norm_tmp$bbp_norm
profiles$sigma_norm <- norm_tmp$sigma_norm
PAR$par_norm <- norm_PAR$par_norm
doxy$doxy_norm <- norm_doxy$doxy_norm

rm(norm_tmp, norm_doxy, norm_PAR)

##########################################
# LET'S PLOT SOMETHING INTERESTING FIRST #
##########################################

# Remove the GS pool from the fit dataframe
fit_no_GS <- fit[!fit$shape %in% 'GS',]
#ggplot(fit_no_GS, aes(x = month, fill = shape)) + geom_bar() + labs(title = 'FIT without GS')

# Only keep GS, GE and G profiles and check which chla_dcm < 1.1*chla_surf
fitggg <- fit[fit$shape %in% c('G','GE','GS'),]
fitggg <- fitggg[which(fitggg$chla_dcm <= 1.1 * fitggg$chla_surf),]
#fit_no_GS <- fit_no_GS[-which(fit_no_GS$juld %in% unique(fitgeg$juld)),]
# Quick histogram to explore 
#ggplot(fitggg, aes(x = month, fill = shape)) + geom_bar() + labs(title = 'G_GE_GS_INF_110%')

# Only keep GS, GE and G profiles and check which chla_dcm < 1.1*chla_surf
fitggg2 <- fit[fit$shape %in% c('G','GE','GS'),]
fitggg2 <- fitggg2[which(fitggg2$chla_dcm >= 1.1 * fitggg2$chla_surf),]
#fit_no_GS <- fit_no_GS[-which(fit_no_GS$juld %in% unique(fitgeg$juld)),]
# Quick histogram to explore 
#ggplot(fitggg2, aes(x = month, fill = shape)) + geom_bar() + labs(title = 'G_GE_GS_SUP_110%')

# # Let's do the same with the GS pool -> see which criteria would be OK for the exclusion threshold
# ggplot(fit, aes(x = month, fill = shape)) + geom_bar()
# fit_thresh <- fit[which(fit$chla_dcm >= 1.1 * fit$chla_surf),]
# ggplot(fit_thresh, aes(x = month, fill = shape)) + geom_bar()
# 
# #sigma ratio
# ggplot(sigma_ratio, aes(x = ratio, fill = factor(month))) + geom_histogram()

# #plot bbp fit
# # same criteria of rejection for GS
# fitbbptest <- fit_BBP[which(fit_BBP$bbp_dbm >= 1.1 * fit_BBP$bbp_surf),]
# ggplot(fit_BBP, aes(x = month, fill = shape)) + geom_bar()
# 
# #correlation
# library(ggpubr)
# fit_thresh$dbm <- fit_BBP$Zmax
# ggscatter(fit_thresh, x = "Zmax", y = "dbm", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")
# #only night
# fit_night <- fit_thresh[fit_thresh$day_night == 'night',]
# ggscatter(fit_night, x = "Zmax", y = "dbm", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")

#####################
# SENSITIVITY STUDY #
#####################

#this dataframe will be used for the creation of the NetCDF files
save_profiles <- profiles

# Let's keep only profiles fitted as "G", "GE" or "GS"
profiles <- clean_remove(profiles, unique(fit[fit$shape %in% c("G", "GE", "GS"),]$juld))

# clean data
fit <- clean_remove(fit, unique(profiles$juld))
density_profiles <- clean_remove(density_profiles, unique(fit$juld))
MLDdf <- clean_remove(MLDdf, unique(fit$juld))
doxy <- clean_remove(doxy, unique(fit$juld))
PAR <- clean_remove(PAR, unique(fit$juld))

# Compute the SIGMA@MAX MLD for each year and each platform
MLDdf_info <- ddply(init_MLD, ~year~platform, summarize,
                    max = max(MLD, na.rm = T),
                    sigmaMaxMLD = sigmaMaxMLD[which.max(MLD)])

#save.image(file="before_sensitivity.RData")
#load(file="before_sensitivity.RData")

# sensitivity
thresholds <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)

sens_ratio <- ldply(as.list(thresholds), function(i){
  
  JULD <- fit[which(fit$chla_dcm >= i * fit$chla_surf),]$juld

  # Compute the sigma_ratio, i.e. SIGMA@DCM/SIGMA@MAX MLD
  sigma_ratioTEST <- ldply(as.list(JULD), function(j){
    tmp <- MLDdf[MLDdf$juld == j,]
    tmp2 <- fit[fit$juld == j,]
    tmp3 <- profiles[profiles$juld == j,]
    sigmaMaxMLD <- MLDdf_info$sigmaMaxMLD[which(MLDdf_info$year == tmp$year & MLDdf_info$platform == tmp$platform)]
    sigmaDCM <- approx(tmp3$depth, tmp3$sigma, tmp2$Zmax)$y #interpolation
    ratio = sigmaDCM/sigmaMaxMLD
    data.frame(ratio = ratio, juld = tmp$juld[1], month = tmp$month[1])
  })
  
  mean_ratio <- mean(sigma_ratioTEST$ratio, na.rm =T)
  sd_ratio <- sd(sigma_ratioTEST$ratio, na.rm =T) 
  var_ratio <- var(sigma_ratioTEST$ratio, na.rm =T)
  data.frame(threshold = i, mean = mean_ratio, std = sd_ratio, var = var_ratio, n = length(JULD))
})



################
# DATA SORTING #
################

# # Let's keep only profiles fitted as "G", "GE" or "GS"
# profiles <- clean_remove(profiles, unique(fit[fit$shape %in% c("G", "GE", "GS"),]$juld))
# fit <- clean_remove(fit, unique(profiles$juld))

# Additional criteria for elimination of false DCMs
# profiles <- clean_remove(profiles, fit[which(fit$chla_dcm > 1.1* fit$chla_surf & fit$Zmax > fit$MLD & fit$month %in% c(3:10)),]$juld) 
profiles <- clean_remove(profiles, fit[which(fit$chla_dcm > 1.1 * fit$chla_surf),]$juld) # let's start with a low threshold (then go to 1.2 or 1.25)
# if we keep the data where Zmax < MLD => check for the ratio at the end.. => if it's good, then we keep it, if it's not good then we have a reason to get rid of them

fit <- clean_remove(fit, unique(profiles$juld))
density_profiles <- clean_remove(density_profiles, unique(fit$juld))
MLDdf <- clean_remove(MLDdf, unique(fit$juld))
doxy <- clean_remove(doxy, unique(fit$juld))
PAR <- clean_remove(PAR, unique(fit$juld))


sigma_ratio <- ldply(as.list(1:length(MLDdf$id)), function(i){
  tmp <- MLDdf[MLDdf$id == i,]
  tmp2 <- fit[fit$id == i,]
  tmp3 <- profiles[profiles$id == i,]
  sigmaMaxMLD <- MLDdf_info$sigmaMaxMLD[which(MLDdf_info$year == tmp$year & MLDdf_info$platform == tmp$platform)]
  sigmaDCM <- approx(tmp3$depth, tmp3$sigma, tmp2$Zmax)$y #interpolation
  ratio = sigmaDCM/sigmaMaxMLD
  data.frame(ratio = ratio, juld = tmp$juld[1], month = tmp$month[1], sigmaMaxMLD = sigmaMaxMLD, sigmaDCM = sigmaDCM)
})

