# PLot trajectories of data?
library(dplyr)
library(plotly)
#######################################
######## BIG FIGURE FOR LIGHT #########
#######################################

#Let's do a function to find the best temporal resolution for our plots

isorho <- function(temp_res, temp_res_iso){ # temp_res = temporal resolution (iso for isolumes, decrease the temporal resolution for those)
  
  time <- seq(from = 1, to = 366, by = temp_res) 
  time_iso <- seq(from = 1, to = 366, by = temp_res_iso) 
  
  length_boxes <- seq(from = 9, to = 15.5, by = 0.2) #15.5 because below that value, we don't care
  nb_boxes <- length(length_boxes)
  
  chla_content <- vector() #for integrating vertical CHLA content (number of elements = number of time steps)
  
  boxesdf <- ldply(as.list(time[1:length(time)-1]), function(j){  # -1 by contruction (i.e. j + temp_res which is already defined)
   
    #check time
    print(j)
    
    #init vectors/arrays
    mean_data <- vector() #mean data for each "box" in the grid ===> FOR MEDIANE  
    number_profile <- vector() #number of profile in each box
    number_obs <- vector()
    #standard_error <- vector()
    chla_pourcentage <- vector()
    q25 <- vector()
    q50 <- vector()
    q75 <- vector()
    chla_content <- vector()
    chlaBBPratio <- vector()
    
    data <- profiles[profiles$DOY >= j & profiles$DOY < j + temp_res,]
    #based on the BOX FRAME ("definition"), we have to cut all data where sigma > 15.5 (15.4) 
    data <- data[-which(data$sigma > 15.4),]
    #remove data without sigma
    data <- data[!is.na(data$sigma),]
    
    # ADD CHLA_SIGMA RATIO
    data$chlaBBPratio <- data$fluo/data$bbp # FLUO must be on top when it's 0. We made the assumptions that BBP is never 0? ==> yes (range(profiles) = [0.00025; 0.04731])
    
    # TO BE DONE: ADD A THRESHOLD FOR THIS SIGMA_RATIO
    
    # integrated content for that period, between min sigma and sigma = 15.4
    #chla_content <- append(chla_content,sum(data$fluo, na.rm = T))
    integrated_content <- sum(data$fluo, na.rm=T)

    for (i in 1:length(length_boxes)-1){
        
        #which data to process?
        index <- which(data$sigma >= length_boxes[i] & data$sigma < length_boxes[i+1]) #take the inf bound but not the sup for boxes separation
        
        # compute data
        mean_i <- mean(data$fluo[index], na.rm = T)
        chla_pc_i <- sum(data$fluo[index], na.rm = T)/integrated_content *100 #normalization
        sd_i <- sd(data$fluo[index])
        N_i <- length(index) #number of observations
        #standard_error_i <- sd_i/sqrt(N_i)
        mean_ratio_chla_bbp_i <- mean(data$chlaBBPratio[index], na.rm = T)
        
        # fill in data arrays
        mean_data[i] <- mean_i
        #standard_error[i] <- standard_error_i
        number_profile[i] <- length(unique(data$id[index]))
        number_obs[i] <- N_i
        chla_content[i] <- sum(data$fluo[index])
        chla_pourcentage[i] <- chla_pc_i
        chlaBBPratio[i] <- mean_ratio_chla_bbp_i
        
        #order fluo data for QUARTILES
        fluo_ordered <- sort(data$fluo[index])
        quart <- quantile(fluo_ordered, c(0.25, 0.5, 0.75), type = 2) #9 methods available, chosen the second one because Q50 = mediation function
        q25[i] <- as.numeric(quart[1])
        q50[i] <- as.numeric(quart[2])
        q75[i] <- as.numeric(quart[3])
    }
    
    # TO CHECK : no data below 15.2 (15.4 has no data i think due to the loop constructions)
    data.frame(sigma_bound = length_boxes[1:nb_boxes-1], mean_data = mean_data, time = rep(j, nb_boxes-1), number_profile = number_profile,
               number_obs = number_obs, q25 = q25, q50 = q50, q75 = q75, chla_pc = chla_pourcentage, chla_content = chla_content, ratioChlaBBP = chlaBBPratio) 
  })
  
  #remove data
  # boxesdf <- boxesdf[-which(is.na(boxesdf$chla_pc)),]
  # boxesdf <- boxesdf[-which(boxesdf$chla_pc == 0),]
  
  boxesdf$number_obs[which(boxesdf$number_obs == 0)] <- NA
  boxesdf$number_profile[which(boxesdf$number_profile == 0)] <- NA
  
  #test plot
  plot_mean_chla <- plot_ly(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$mean_data, nrow = nb_boxes-1, ncol = length(time)), type = "contour",
               contours = list(start = 0, end = 1.5, size = 0.2))
  
  plot_chla_pourcentage <- plot_ly(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$chla_pc, nrow = nb_boxes-1, ncol = length(time)), type = "contour",
               contours = list(start = 0, end = 30, size = 2))
  
  plot_CHLA_BBP_ratio <- plot_ly(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$ratioChlaBBP, nrow = nb_boxes-1, ncol = length(time)), type = "contour",
                                   contours = list(start = 0, end = 600, size = 50))
  
  # boxesdf$abundance <- boxesdf$chla_pc * boxesdf$ratioChlaBBP # test
  # plot_ly(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$abundance, nrow = nb_boxes-1, ncol = length(time)), type = "contour",
  #         contours = list(start = 0, end = 10000, size = 500))
  
  #add light curves PAR_0, take all data available
  par_density <- profiles[profiles$juld %in% PAR$juld,]
  par_density <- transform(par_density, id=as.numeric(factor(juld)))
  
  PAR2 <- PAR[PAR$juld %in% par_density$juld,]
  PAR2 <- transform(PAR2, id=as.numeric(factor(juld)))
  
  parPC <- ldply(as.list(1:length(unique(PAR2$id))), function(i){ #PC for PourCent
    
    #code pourri ici mais ça marche et j'ai pas le temps..
    
    tmp <- PAR2[PAR2$id == i,]
    par0 <- tmp$par[1]
    pc <- vector()
    for (i in 1:length(tmp$depth)){
      pc[i] <- tmp$par[i]/par0
    }
    
    data.frame(parpc = pc)
    
  })
  
  #put back in par_dcm the par in %
  PAR2 <- PAR2[order(PAR2$id),]
  PAR2$parpc <- parPC$parpc
  
  #interpolate sigma on parPC
  a <- ldply(as.list(1:length(unique(par_density$id))), function(i){
    
    tmp <- par_density[par_density$id == i,]
    tmp2 <- PAR2[PAR2$id == i,]
    par_density2 <-approx(tmp$depth, tmp$sigma, tmp2$depth)$y #interpoler bcp plus simple que le truc que j'ai pu faire avant bordel
    tmp2 <- cbind(tmp2, par_density2)
    
    data.frame(depth=tmp2$depth, juld=tmp2$juld, par = tmp2$par, par_sigma = par_density2)
    
    #i <- i + 1
  })
  
  PAR2$sigma_par <- a$par_sigma
  
  #add light curves #approx, not interpolated
  b10 <- ldply(as.list(1:length(unique(par_density$id))), function(i){
    
    tmp <- par_density[par_density$id == i,]
    curve <- 0.1 # courbe des 10% du PAR0
    tmp2 <- PAR2[PAR2$id == i,]
    ind <- which.min(abs(curve-tmp2$parpc))
    sigma_pc <- tmp$sigma[ind] 
    
    data.frame(sigma_pc = sigma_pc, DOY = tmp$DOY[1])
    
  })
  
  b1 <- ldply(as.list(1:length(unique(par_density$id))), function(i){
    
    tmp <- par_density[par_density$id == i,]
    tmp2 <- PAR2[PAR2$id == i,]
    curve <- 0.01 # courbe des 1% du PAR0
    ind <- which.min(abs(curve-tmp2$parpc))
    sigma_pc <- tmp$sigma[ind] 
    
    data.frame(sigma_pc = sigma_pc, DOY = tmp$DOY[1])
  })
  
  mean_iso_10pc <- vector() #mean data for each "box" in the grid
  std_iso_10pc <- vector() #number of profile in each box
  mean_iso_1pc <- vector()
  std_iso_1pc <- vector()
  q25iso1pc <- vector()
  q50iso1pc <- vector()
  q75iso1pc <- vector()
  q25iso10pc <- vector()
  q50iso10pc <- vector()
  q75iso10pc <- vector()
  
  isolumes_times <- ldply(as.list(time_iso[1:length(time_iso)-1]), function(j){  #3:10 because March -> October for DCMs (le -1 est là pour ne pas s'emmerder avec les détails de fin de tableau)
    #init
    print(j)
    data_iso_1pc <- b1[b1$DOY >= j & b1$DOY < j + temp_res_iso,]
    data_iso_10pc <- b10[b10$DOY >= j & b1$DOY < j + temp_res_iso,]
    
    mean_iso_1pc <- append(mean_iso_1pc, mean(data_iso_1pc$sigma_pc, na.rm = T))
    std_iso_1pc <- append(std_iso_1pc, sd(data_iso_1pc$sigma_pc, na.rm = T))
    
    mean_iso_10pc <- append(mean_iso_10pc, mean(data_iso_10pc$sigma_pc, na.rm = T))
    std_iso_10pc <- append(std_iso_10pc, sd(data_iso_10pc$sigma_pc, na.rm = T))
    
    #order before computing quartiles
    #remove NA's
    #data_iso_10pc <- data_iso_10pc[-which(is.na(data_iso_10pc$sigma_pc)),]
    #data_iso_1pc <- data_iso_1pc[-which(is.na(data_iso_1pc$sigma_pc)),]
    orderiso10pc <- sort(data_iso_10pc$sigma_pc, na.last=NA)
    orderiso1pc <- sort(data_iso_1pc$sigma_pc, na.last=NA)
    
    #quartiles
    #10pc
    q25iso10pc <- append(q25iso10pc, as.numeric(quantile(orderiso10pc, 0.25, type = 2)))
    q50iso10pc <- append(q50iso10pc, as.numeric(quantile(orderiso10pc, 0.50, type = 2))) #mediane
    q75iso10pc <- append(q75iso10pc, as.numeric(quantile(orderiso10pc, 0.75, type = 2)))
    #1pc
    q25iso1pc <- append(q25iso1pc, as.numeric(quantile(orderiso1pc, 0.25, type = 2)))
    q50iso1pc <- append(q50iso1pc, as.numeric(quantile(orderiso1pc, 0.50, type = 2)))
    q75iso1pc <- append(q75iso1pc, as.numeric(quantile(orderiso1pc, 0.75, type = 2)))
    
    data.frame(time = j, mean_iso_1pc = mean_iso_1pc, mean_iso_10pc = mean_iso_10pc, q25iso10pc = q25iso10pc, q50iso10pc = q50iso10pc, q75iso10pc = q75iso10pc,
               q25iso1pc = q25iso1pc, q50iso1pc = q50iso1pc, q75iso1pc = q75iso1pc)
  })
  
  TOTALPLOT <- plot_ly() %>%
    #add_trace(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$mean_data, nrow = nb_boxes-1, ncol = length(time)),
    #          type = 'contour', contours = list(start = 0, end = 1.5, size = 0.2)) %>%
    add_trace(x = seq(from = 1, to = 366, by = temp_res), y = -length_boxes[1:length(length_boxes-1)], z = matrix(boxesdf$chla_pc, nrow = nb_boxes-1, ncol = length(time)),
              type = 'contour', contours = list(start = 0, end = 30, size = 2)) %>%
    add_lines(x = time, y = -13.5, color = I("darkgreen"), name = "sigma_limit_13.5") %>%
    add_lines(x = isolumes_times$time, y = -isolumes_times$q25iso10pc, type ='scatter', name = "iso 10% Q25", line = list(color = "black", dash = "dash")) %>% #, error_y = ~list(array = isolumes_times$std_iso_10pc, color = "black")) %>%
    add_lines(x = isolumes_times$time, y = -isolumes_times$q50iso10pc, type ='scatter', name = "iso 10% Q50", line = list(color = "black")) %>% #, error_y = ~list(array = isolumes_times$std_iso_1pc, color = "red"))
    add_lines(x = isolumes_times$time, y = -isolumes_times$q75iso10pc, type ='scatter', name = "iso 10% Q75", line = list(color = "black", dash = "dot")) %>%#, error_y = ~list(array = isolumes_times$std_iso_1pc, color = "red"))
    add_lines(x = isolumes_times$time, y = -isolumes_times$q25iso1pc, type ='scatter', name = "iso 1% Q25", line = list(color = "red", dash = "dash")) %>%#, error_y = ~list(array = isolumes_times$std_iso_1pc, color = "red"))
    add_lines(x = isolumes_times$time, y = -isolumes_times$q50iso1pc, type ='scatter', name = "iso 1% Q50", line = list(color = "red")) %>%#, error_y = ~list(array = isolumes_times$std_iso_1pc, color = "red"))
    add_lines(x = isolumes_times$time, y = -isolumes_times$q75iso1pc, type ='scatter', name = "iso 1% Q75", line = list(color = "red", dash = "dot")) #, error_y = ~list(array = isolumes_times$std_iso_1pc, color = "red"))

  
  return(TOTALPLOT)
}

# Separation of the integrated content ===> roughly done
limits <- 13.5 # Visually, maybe we can put 13
test <- boxesdf
test_up <- boxesdf[boxesdf$sigma_bound <= limits,]
test_down <- boxesdf[boxesdf$sigma_bound > limits,]

#check between DOY = 150 and 240 (~ Summer)
library(data.table)
total_chla <- sum(boxesdf[boxesdf$time %between% c(150,250),]$chla_content)
total_chla_summer_up <- sum(test_up[test_up$time %between% c(150,250),]$chla_content)
total_chla_summer_down <- sum(test_down[test_down$time %between% c(150,250),]$chla_content)

# Integrated content over time (below and above the sigma limit [13.5)])

integrated_content_up <- vector()
integrated_content_down <- vector()
all_content <- vector()
all_type <- vector()
all_time <- vector()

integrated_chla_time <- ldply(as.list(time[1:length(time)-1]), function(j){ 
  
  print(j)
  data_up <- test_up[test_up$time >= j & test_up$time < j + temp_res,]
  data_down <- test_down[test_down$time >= j & test_down$time < j + temp_res,]
  data <- test[test$time >= j & test$time < j + temp_res,]
  
  integrated_content_up <- append(integrated_content_up, sum(data_up$chla_content, na.rm=T)/sum(data$chla_content, na.rm=T)*100)
  integrated_content_down <- append(integrated_content_down, sum(data_down$chla_content, na.rm=T)/sum(data$chla_content, na.rm=T)*100)
  up <- "UP"
  down <- "DOWN"
  
  all_content <- append(integrated_content_up, integrated_content_down)
  # up <- rep("UP", length(time))
  # down <- rep("DOWN", length(time))
  all_type <- append(up, down)
  all_time <- append(j, j)
  
  data.frame(time = all_time, content = all_content, type = all_type)
  })

down <- integrated_chla_time[integrated_chla_time$type == "DOWN",]
up <- integrated_chla_time[integrated_chla_time$type == "UP",]

ggplot(down, aes(x = time, y = content)) + geom_area(fill ="red", alpha = .5) +
  geom_area(data = up, fill = "green", alpha=.5)
    