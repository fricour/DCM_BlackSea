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
  data.frame(chla_norm = chla_norm, bbp_norm = bbp_norm)
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

# insert nomalized CHLA and BBP into PROFILES
profiles$fluo_norm <- norm_tmp$chla_norm
profiles$bbp_norm <- norm_tmp$bbp_norm
PAR$par_norm <- norm_PAR$par_norm
doxy$doxy_norm <- norm_doxy$doxy_norm

rm(norm_tmp, norm_doxy, norm_PAR)

#TBD : what is following, the application and the thing asked by Arthur

# MLDdf3 <- ddply(MLDdf_alltime, ~year~Platform, summarize,
#                 max = max(MLD, na.rm = T),
#                 sigmaMaxMLD = sigmaMaxMLD[which.max(MLD)])
# 
# argodf <- cbind(argodf, sigma_dcm_from_depth_dcm)
# colnames(argodf)[20] <- "sigma_dcm"
# 
# sigma_ratio <- ldply(as.list(1:length(MLDdf$id)), function(i){
#   tmp <- MLDdf[MLDdf$id == i,]
#   tmp2 <- argodf[argodf$id == i,]
#   sigmaMaxMLD <- MLDdf3$sigmaMaxMLD[which(MLDdf3$year == tmp$year & MLDdf3$Platform == tmp$Platform)]
#   ratio = tmp2$sigma_dcm/sigmaMaxMLD
#   data.frame(ratio = ratio, juld = tmp$juld[1])
# })

# 
# ggplot(sigma_ratio, aes(ratio)) +  geom_histogram(data = sigma_ratio, col = "black")+
#   theme(legend.position="none") +
#   ylab("Number of profiles") + xlab(expression(paste(sigma["DCM"],"/",sigma["MLD-MAX"]))) +
#   theme(text = element_text(size=13))
