#Function to simulate a single spatial field - returns df of coordinates and values

sim_field <- function(mask, #vector of polygons of countries/boundaries
                      vgm_model = "Gau", #variogram model type - "Exp", "Sph", or "Gau"
                      nugget = 0, #variance at 0 distance
                      sill = 1, #asymptotic max variance
                      range = 150000, #distance at which 95% of max variance reached
                      k = 4, #neighbors to consider when simulating 
                      n_sample = 1000 #number of points to sample for simulation
                      
){
  
  ###
  # Prepare geodata
  ###
  
  #sample points from the vector
  pointsample <- terra::spatSample(mask, size = n_sample, method = "regular")
  
  #get coordinates
  xy <- as.data.frame(crds(pointsample, df = TRUE))
  
  ###
  # Define spatial model and predict
  ###
  
  spat_model <- gstat(formula = z ~ 1, #single variable no dependence
                      locations = ~x+y,
                      dummy = TRUE, #unconditional simulation
                      beta = nugget, 
                      model = vgm(model = vgm_model,
                                  psill = sill, 
                                  range = range
                      ),
                      nmax = k
  )
  
  spat_sim <- predict(spat_model, 
                      newdata = xy, 
                      nsim = 1,
                      debug.level = 0) |> 
    mutate(t1 = sim1) |> 
    select(-sim1)
  
  return(spat_sim)
}
