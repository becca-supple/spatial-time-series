###
# Load data
###

GBFS_complete <- read.csv(file = "GBFS/GBFS_11_36.csv")
load(file = "GBFS/nb_complete.RData")

#set significance levels and weights
alphas <- c(0.001, 0.01, 0.025, 0.05, 0.1)
omega <- c(100/186, 50/186, 25/186, 10/186, 1/186)

###
# Fit by species, by alpha value:
###

#House Sparrow
GBFS_housesparrow <- GBFS_complete |> 
  select(county, year, subrur, temp, frost, sparr, housp)

housp_ind <- list()
housp_oriented <- list()

for(a in seq(1, 5)){
  
  housp_ind[[a]] <- conditional_independencies(data = GBFS_housesparrow,
                                               max_lag = 2,
                                               factor = "county",
                                               time = "year",
                                               family = c("binomial", "gaussian", "gaussian", "tw", "tw"),
                                               neighborhood = nb_complete,
                                               alpha = alphas[a])
  
  housp_oriented[[a]] <- orient_edges(housp_ind[[a]])
  
}

housp_swag <- edge_dist_weight(housp_oriented, weights = omega)

#Starling

GBFS_starling <- GBFS_complete |> 
  select(county, year, subrur, temp, frost, sparr, starl)

starl_ind <- list()
starl_oriented <- list()

for(a in seq(1, 5)){

  print(Sys.time())
  print(alphas[a])
  
  starl_ind[[a]] <- conditional_independencies(data = GBFS_starling,
                                               max_lag = 2,
                                               factor = "county",
                                               time = "year",
                                               family = c("binomial", "gaussian", "gaussian", "tw", "tw"),
                                               neighborhood = nb_complete,
                                               alpha = alphas[a])
  
  starl_oriented[[a]] <- orient_edges(starl_ind[[a]])
  
}

starl_swag <- edge_dist_weight(starl_oriented, weights = omega)

#Blue Tit

GBFS_bluetit <- GBFS_complete |> 
  select(county, year, subrur, temp, frost, sparr, bluti)

bluti_ind <- list()
bluti_oriented <- list()

for(a in seq(1, 5)){

  print(Sys.time())
  print(alphas[a])
  
  bluti_ind[[a]] <- conditional_independencies(data = GBFS_bluetit,
                                               max_lag = 2,
                                               factor = "county",
                                               time = "year",
                                               family = c("binomial", "gaussian", "gaussian", "tw", "tw"),
                                               neighborhood = nb_complete,
                                               alpha = alphas[a])
  
  bluti_oriented[[a]] <- orient_edges(bluti_ind[[a]])
  
}

bluti_swag <- edge_dist_weight(bluti_oriented, weights = omega)

#Greenfinch

GBFS_greenfinch <- GBFS_complete |> 
  select(county, year, subrur, temp, frost, sparr, grefi)

grefi_ind <- list()
grefi_oriented <- list()

for(a in seq(1, 5)){

  print(Sys.time())
  print(alphas[a])
  
  grefi_ind[[a]] <- conditional_independencies(data = GBFS_greenfinch,
                                               max_lag = 2,
                                               factor = "county",
                                               time = "year",
                                               family = c("binomial", "gaussian", "gaussian", "tw", "tw"),
                                               neighborhood = nb_complete,
                                               alpha = alphas[a])
  
  grefi_oriented[[a]] <- orient_edges(grefi_ind[[a]])
  
}

grefi_swag <- edge_dist_weight(grefi_oriented, weights = omega)


