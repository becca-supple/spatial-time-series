
sim_system_linear <- function(mask, #vector of polygons of countries/boundaries
                              area_ID, #name of object that stores the area IDs in the mask
                              t_steps, #number of time steps
                              autocorrelation = 0.1, #autocorrelation between years
                              vgm_data = "Gau", #variogram model type  for data - "Exp", "Sph", or "Gau"
                              nugget = 0, #variance at 0 distance
                              sill = 1, #asymptotic max variance
                              range = 150000, #distance at which 95% of max variance reached
                              k = 4, #neighbors to consider when simulating 
                              n_sample = 2500, #number of points to sample for simulation
                              strength = c("none", "low", "med", "high"), #strength of relationships
                              vgm_rel = "Gau", #variogram model type for relationship parameters
                              nugget_rel = 0, #parameter variance at 0 distance
                              sill_rel = 0.001, #asymptotic max variance for parameters of relationships
                              range_rel = 4000000, #distance at which 95% of max parameter variance reached
                              noise = 0.1, #random noise
                              return_rast = FALSE #return raster data
){
  
  ###
  # Simulate Agriculture
  ###
  
  #Set up blank list
  agr <- list()
  
  #Initialize spatial field
  agr_field <- sim_field(mask = mask,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill+5,
                         range = range,
                         k = k,
                         n_sample = n_sample)
  
  n_locations <- nrow(agr_field)
  xy <- agr_field[, c(1, 2)]
  
  #Initialize time series
  agr[[1]] <- agr_field[["t1"]]
  
  for(t in seq(2, t_steps + 1)){
    
    agr[[t]] <- numeric(n_locations)
    
    epsilon <- sim_field(mask = mask,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill,
                         range = range,
                         k = k,
                         n_sample = n_sample)[["t1"]]
    
    agr[[t]] <- autocorrelation * agr[[t-1]] + sqrt(1 - autocorrelation^2) * epsilon
    
  }
  
  #Transform data
  agr <- lapply(agr, function(i){
    i + rnorm(length(i), 0, noise)
  })
  agr <- lapply(agr, pnorm) #probit link to make data in [0,1]
  
  
  ###
  # Simulate Fencing
  ###
  
  #Set up blank list
  fence <- list()
  
  #Spatial variance for fence-agr relationships
  beta_fence <- sim_field(mask = mask,
                          vgm_model = vgm_rel,
                          nugget = nugget_rel,
                          sill = sill_rel,
                          range = range_rel,
                          k = k,
                          n_sample = n_sample)
  
  #Vector of parameters
  beta <- ifelse(strength == "none", 0, 
                 ifelse(strength == "low", 0.5, 
                        ifelse(strength == "med", 1, 2)))
  
  betas <- beta + beta_fence$t1
  
  #Initialize time series
  epsilon <- sim_field(mask = mask,
                       vgm_model = vgm_data,
                       nugget = nugget,
                       sill = sill,
                       range = range,
                       k = k,
                       n_sample = n_sample)[["t1"]]
  
  fence[[1]] <- beta * agr[[1]] + 
    sqrt(1 - autocorrelation^2) * epsilon + 
    rnorm(n_locations, 0, noise)
  
  for(t in seq(2, t_steps + 1)){
    
    fence[[t]] <- numeric(n_locations)
    
    epsilon <- sim_field(mask = mask,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill,
                         range = range,
                         k = k,
                         n_sample = n_sample)[["t1"]]
    
    fence[[t]] <- autocorrelation * fence[[t-1]] + 
      betas * agr[[t]] + 
      sqrt(1 - autocorrelation^2) * epsilon
    
  }
  
  #center and scale data
  fence <- lapply(fence, function(i){
    i + rnorm(length(i), 0, noise)
  })
  fence <- lapply(fence, function(i){
    (i - mean(i))/(sd(i))
  })
  
  ###
  # Simulate Fertilizer
  ###
  
  #Set up blank list
  fert <- list()
  
  #Spatial variance for fert-agr relationships
  beta_fert <- sim_field(mask = mask,
                         vgm_model = vgm_rel,
                         nugget = nugget_rel,
                         sill = sill_rel,
                         range = range_rel,
                         k = k,
                         n_sample = n_sample)
  
  #Vector of parameters
  betas <- beta + beta_fert$t1
  
  #Initialize time series
  epsilon <- sim_field(mask = mask,
                       vgm_model = vgm_data,
                       nugget = nugget,
                       sill = sill,
                       range = range,
                       k = k,
                       n_sample = n_sample)[["t1"]]
  
  fert[[1]] <- beta * agr[[1]] + sqrt(1 - autocorrelation^2) * epsilon + rnorm(n_locations, 0, noise)

  
  for(t in seq(2, t_steps + 1)){
    
    fert[[t]] <- numeric(n_locations)
    
    epsilon <- sim_field(mask = mask,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill,
                         range = range,
                         k = k,
                         n_sample = n_sample)[["t1"]]
    
    fert[[t]] <- autocorrelation * fert[[t-1]] + 
      betas * agr[[t]] + 
      sqrt(1 - autocorrelation^2) * epsilon
    
  }
  
  #center and scale data
  fert <- lapply(fert, function(i){
    i + rnorm(length(i), 0, noise)
  })
  fert <- lapply(fert, function(i){
    (i - mean(i))/(sd(i))
  })
  
  ###
  # Simulate Occupancy
  ###
  
  #Set up blank list
  occ <- list()
  p_occ <- list()
  
  #Spatial variance for fert-occ relationships
  beta_fert_occ <- sim_field(mask = mask,
                             vgm_model = vgm_rel,
                             nugget = nugget_rel,
                             sill = sill_rel,
                             range = range_rel,
                             k = k,
                             n_sample = n_sample)
  
  #Vector of parameters
  betas_fert <- beta/2 + beta_fert_occ$t1
  
  #Spatial variance for agr-occ relationships
  beta_agr_occ <- sim_field(mask = mask,
                            vgm_model = vgm_rel,
                            nugget = nugget_rel,
                            sill = sill_rel,
                            range = range_rel,
                            k = k,
                            n_sample = n_sample)
  
  #Vector of parameters
  betas_agr <- beta + beta_agr_occ$t1
  
  #Initialize field
  occ[[1]] <- sim_field(mask = mask,
                        vgm_model = vgm_data,
                        nugget = nugget,
                        sill = sill,
                        range = range,
                        k = k,
                        n_sample = n_sample)[["t1"]]
  
  occ[[1]] <- sapply(occ[[1]], function(i){
    i + rnorm(1, 0, noise) 
  })
  
  #Sim forward
  for(t in seq(2, t_steps + 1)){

    occ[[t]] <- numeric(n_locations)
    
    epsilon <- sim_field(mask = mask,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill,
                         range = range,
                         k = k,
                         n_sample = n_sample)[["t1"]]
    
    occ[[t]] <- autocorrelation * occ[[t-1]] + 
      betas_agr * agr[[t]] + 
      (-1) * betas_fert * fert[[t-1]] +
      sqrt(1 - autocorrelation^2) * epsilon
    
  }
  
  #transform
  occ <- lapply(occ, function(i){
    i + rnorm(length(i), 0, noise)
  })
  occ <- lapply(occ, pnorm) 
  
  ###
  # Summarize by Area
  ###
  
  if(return_rast){
    rast_stack <- list()
  }
  
  data_short <- data.frame(area_ID = unique(mask[[area_ID]]))
  
  
  ##Summarize agriculture##
  
  for(t in seq(2, t_steps + 1)){
    
    #attach data to coordinates
    spat_sim <- cbind(xy, agr[[t]])
    
    #make raster
    rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
    rast_sim <- mask(rast_sim, mask)
    
    if(return_rast){
      rast_stack[[t-1]] <- c(rast_sim)
    }
    
    #extract mean by country
    country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)
    
    data_short[[paste0("T", t - 1)]] <- country_means[,2]
    
  }
  
  #Pivot data
  data_long <- pivot_longer(data_short,
                            cols = paste0("T", seq(1, t_steps)),
                            names_to = "Time",
                            values_to = "Agr")
  
  ##Summarize fencing##
  
  for(t in seq(2, t_steps + 1)){
    
    #attach data to coordinates
    spat_sim <- cbind(xy, fence[[t]])
    
    #make raster
    rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
    rast_sim <- mask(rast_sim, mask)
    
    if(return_rast){
      rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
    }
    
    #extract mean by country
    country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)
    
    data_short[[paste0("T", t - 1)]] <- country_means[,2]
    
  }
  
  #Pivot data
  fence_long <- pivot_longer(data_short,
                             cols = paste0("T", seq(1, t_steps)),
                             names_to = "Time",
                             values_to = "Fence")
  
  data_long[["Fence"]] <- fence_long[["Fence"]]
  
  ##Summarize fertilizer##
  
  for(t in seq(2, t_steps + 1)){
    
    #attach data to coordinates
    spat_sim <- cbind(xy, fert[[t]])
    
    #make raster
    rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
    rast_sim <- mask(rast_sim, mask)
    
    if(return_rast){
      rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
    }
    
    #extract mean by country
    country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)
    
    data_short[[paste0("T", t - 1)]] <- country_means[,2]
    
  }
  
  #Pivot data
  fert_long <- pivot_longer(data_short,
                            cols = paste0("T", seq(1, t_steps)),
                            names_to = "Time",
                            values_to = "Fert")
  
  data_long[["Fert"]] <- fert_long[["Fert"]]
  
  ##Summarize occupancy##
  
  for(t in seq(2, t_steps + 1)){
    
    #attach data to coordinates
    spat_sim <- cbind(xy, occ[[t]])
    
    #make raster
    rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
    rast_sim <- mask(rast_sim, mask)
    
    if(return_rast){
      rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
      names(rast_stack[[t-1]]) <- c("Agr", "Fence", "Fert", "Occ")
    }
    
    #extract mean by country
    country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)
    
    data_short[[paste0("T", t - 1)]] <- country_means[,2]
    
  }
  
  #Pivot data
  occ_long <- pivot_longer(data_short,
                           cols = paste0("T", seq(1, t_steps)),
                           names_to = "Time",
                           values_to = "Occ")
  
  data_long[["Occ"]] <- occ_long[["Occ"]]
  
  ### 
  # Format and return
  ###
  
  data_long <- data_long |> 
    filter(!is.na(Agr)) |> #remove countries with NA (out of frame)
    mutate(Time = as.numeric(gsub("T", "", Time))) #numeric time steps
  
  if(return_rast){
    list2return <- list(data = data_long, rasters = rast_stack)
    return(list2return)
  }else{return(data_long)}
  
}
