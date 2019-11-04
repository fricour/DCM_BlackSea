#########################
##### BEGIN THE FIT #####
#########################

fit_BBP <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id==i,]
  
  # restrain data to the oxygenated layer (to avoid fitting the second BBP peak at the anoxic interface)
  tmp_doxy <- doxy[doxy$id == i,]
  
  # TEST : % of the oxygen content at the surface? Or > 10 mumol/kg?
  doxy_threshold <- 10
  depth_doxy <- tmp_doxy$depth[which.max(tmp_doxy$doxy < doxy_threshold)]
  if(depth_doxy < 1){ # can happen if doxy > 10 for the whole profile (when it's a shallower profile than usual) # happened only when threshold was fixed at 5
    depth_doxy <- max(tmp$depth)
  }
  tmp <- tmp[tmp$depth <= depth_doxy,]
  
  #data needed for R^2 and R^2 adjusted
  mean_bbp <- mean(tmp$bbp, na.rm=T) #mean bbp over all profile
  ss_tot <- sum((tmp$bbp - mean_bbp)^2, na.rm=T) #for R^2 computation
  n <- length(tmp$depth) #number of data points
  #Arthur's suggestion => Fsurf should be a mean of the fluo in the MLD ==> here then it's the same idea for BBP
  index_MLD <- which.min(tmp$depth <= MLDdf$MLD[i]) 
  
  ##############
  # PARAMETERS #
  ##############
  
  #Parameters estimation
  Fsurf <- mean(tmp$bbp[1:index_MLD], na.rm = T)
  Zdemi <- tmp$depth[which.max(tmp$bbp <=  Fsurf/2)]
  if(MLDdf$MLD[i] > Zdemi){
    Zdemi <- MLDdf$MLD[i]
  }
  maxindex <- which.max(tmp$bbp)
  Zmax <- tmp$depth[maxindex]
  Fmax <- tmp$bbp[maxindex]
  dz <- 5 # First guess, cannot do better
  s <- -0.01 #guess without explication like dz, it is fixed and adapted if initial values do not work
  
  #############################
  # GAUSSIAN-EXPONENTIAL (GE) #
  #############################
  
  res_GE <- tryCatch(nlsLM(bbp ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                           start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                           data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                    minFactor = 1/2048, warnOnly=T)),
                     error=function(e) e)
  
  
  if(inherits(res_GE, "error")){# handle error if initial dz is not good
    dz <- 10
    res_GE <- tryCatch(nlsLM(bbp ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                             start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                             data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                      minFactor = 1/2048, warnOnly=T)),
                       error=function(e) e)
  }
  
  #R^2 and R^2 adjusted
  ge <- fgauss_expo(tmp$depth,coef(res_GE)["Fsurf"], coef(res_GE)["Zdemi"], coef(res_GE)["Fmax"], coef(res_GE)["Zmax"], coef(res_GE)["dz"])
  ss_res_ge <- sum((ge - tmp$bbp)^2, na.rm=T)  
  rcoef_ge <- 1-(ss_res_ge/ss_tot)
  rcoef_ge <- 1 - (1-rcoef_ge)*(n-1)/(n-5-1) # final R^2 ADJUSTED
  
  ###############
  # SIGMOID (S) #
  ###############
  
  res_S <- tryCatch(nlsLM(bbp ~ fsigmoid(depth, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                   minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  
  if(inherits(res_S, "error")){#added due to a bug
    Zdemi <- 30
    #s <- -0.01
    res_S <- tryCatch(nlsLM(bbp ~ fsigmoid(depth, Fsurf, Zdemi, s),
                            start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                            data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                     minFactor = 1/2048, warnOnly=T)),
                      error=function(e) e)
  }
  
  
  sigmoid <- fsigmoid(tmp$depth, coef(res_S)["Fsurf"], coef(res_S)["Zdemi"], coef(res_S)["s"])
  ss_res_sigmoid <- sum((sigmoid - tmp$bbp)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)               
  rcoef_sigmoid <- 1 - (1-rcoef_sigmoid)*(n-1)/(n-3-1)
  
  ################
  # GAUSSIAN (G) #
  ################
  
  res_G <- tryCatch(nlsLM(bbp ~ fgauss(depth, Fmax, Zmax, dz),
                          start = c(Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                   minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  
  gauss <- fgauss(tmp$depth, coef(res_G)["Fmax"], coef(res_G)["Zmax"], coef(res_G)["dz"])
  ss_res_gauss <- sum((gauss - tmp$bbp)^2, na.rm=T)                
  rcoef_gauss <- 1-(ss_res_gauss/ss_tot)  
  rcoef_gauss <- 1 - (1-rcoef_gauss)*(n-1)/(n-3-1)
  
  ###################
  # EXPONENTIAL (E) #
  ###################
  
  res_E <- tryCatch(nlsLM(bbp ~ fexpo(depth, Fsurf, Zdemi),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi),
                          data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                   minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  
  expo <- fexpo(tmp$depth, coef(res_E)["Fsurf"], coef(res_E)["Zdemi"])
  ss_res_expo <- sum((expo - tmp$bbp)^2, na.rm=T)
  rcoef_expo <- 1-(ss_res_expo/ss_tot)
  rcoef_expo <- 1 - (1-rcoef_expo)*(n-1)/(n-2-1)
  
  #########################
  # GAUSSIAN-SIGMOID (GS) #
  #########################
  
  res_GS <- tryCatch(nlsLM(bbp ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                                 Fmax, Zmax, dz),
                           start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                           data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                    minFactor = 1/2048, warnOnly=T)),
                     error=function(e) e)
  
  if(inherits(res_GS, "error")){
    Zdemi <- 30#added due to a bug
    s <- -0.1 #same as before, guess without explainations
    res_GS <- tryCatch(nlsLM(bbp ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                                   Fmax, Zmax, dz),
                             start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                             data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                      minFactor = 1/2048, warnOnly=T)),
                       error=function(e) e)
  }
  
  if(inherits(res_GS, "error")){
    dz <- 10
    s <- -0.1
    res_GS <- tryCatch(nlsLM(bbp ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                                   Fmax, Zmax, dz),
                             start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                             data=tmp[,c(1,4)], control = nls.control(maxiter=1000,
                                                                      minFactor = 1/2048, warnOnly=T)),
                       error=function(e) e)
  }
  
  gs <- fgauss_sigmoid(tmp$depth, coef(res_GS)["Fsurf"], coef(res_GS)["Zdemi"], coef(res_GS)["s"], coef(res_GS)["Fmax"], coef(res_GS)["Zmax"], coef(res_GS)["dz"])
  ss_res_gs <- sum((gs - tmp$bbp)^2, na.rm=T)
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
    bbp_at_dbm <- NA
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
    bbp_at_dbm <- NA
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
    bbp_at_dbm <- NA
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
    bbp_at_dbm <- tmp$bbp[which.min(abs(tmp$depth - Zmax))] 
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
    bbp_at_dbm <- tmp$bbp[which.min(abs(tmp$depth - Zmax))]
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
    bbp_at_dbm <- tmp$bbp[which.min(abs(tmp$depth - Zmax))]
    Fsurf <- coef(res_GS)["Fsurf"]
    Zdemi <- coef(res_GS)["Zdemi"]
    s <- coef(res_GS)["s"]
    Fmax <- coef(res_GS)["Fmax"]
    Zmax <- coef(res_GS)["Zmax"]
    dz <- coef(res_GS)["dz"]
  }
  
  show(i)
  
  data.frame(shape = shape, r2adj = r2adj, Zmax = Zmax, id = tmp$id[1], juld = tmp$juld[1],
             lon = tmp$lon[1], lat = tmp$lat[1],
             month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$platform[1], bbp_surf = tmp$bbp[1], #bbp_surf is the bbp value at the surface of the profile
             day_night = tmp$day_night[1], MLD = MLDdf$MLD[i], bbp_dbm = bbp_at_dbm, Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, dz = dz) 
})
