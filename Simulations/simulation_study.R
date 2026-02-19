###
# Set up parameter dataframe
###

time_length <- c(10, 20, 30)
noise <- c(0.01, 0.25, 0.5, 0.75, 1)
autocorrelation <- c(0.01, 0.25, 0.5, 0.75)
strength <- c("none", "low", "med", "high")
range <- c(100000, 150000, 300000)

params <- expand.grid(time_length, noise, autocorrelation, 
                      strength, range)

colnames(params) <- c("time_length", "noise", 
                      "autocorrelation", "strength", "range")

#how many systems to simulate per scenario
n_reps <- 3

#how many times to allow failure before moving on
max_fails <- 1

###
# Simulate systems
###

test_systems_l <- list()
seed_list <- list()

start <- Sys.time()
for(i in 1:nrow(params)){
  
  tl <- params$time_length[i]
  ns <- params$noise[i]
  ac <- params$autocorrelation[i]
  st <- params$strength[i]
  rg <- params$range[i]
  
  for(rep in seq(1, n_reps)){
    
    #progress
    print(paste0(i, "_", rep))
    
    #set a random seed
    set.seed(NULL)
    seed <- seed_list[[paste0(i, "_", rep)]]
    set.seed(seed)
    
    #simulate
    test_system <- sim_system_linear(mask = eur_boundaries_filtered,
                                     area_ID = "Country",
                                     t_steps = tl,
                                     strength = st,
                                     autocorrelation = ac,
                                     noise = ns,
                                     range = rg,
                                     return_rast = FALSE)
    
    test_systems_l[[paste0(i, "_", rep)]] <- test_system
  }
}

###
# Run causal discovery algorithm
###

test_inds_l <- list()
test_oriented_l <- list()

for(i in 1:nrow(params)){
  
  for(rep in seq(1, n_reps)){
    
    #progress
    print(paste0(i, "_", rep))

    test_inds_l[[paste0(i, "_", rep)]] <- conditional_independencies(data = test_systems_l[[paste0(i, "_", rep)]],
                                                                     max_lag = 2,
                                                                     factor = "Country",
                                                                     time = "Time",
                                                                     family = c("betar", "gaussian", "gaussian", "betar"),
                                                                     neighborhood = nb,
                                                                     modeltype = "separated")
    
    test_oriented_l[[paste0(i, "_", rep)]] <- orient_edges(test_inds_l[[paste0(i, "_", rep)]])
    
  }
}
