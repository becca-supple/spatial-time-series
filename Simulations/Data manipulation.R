#load necessary packages
library(dplyr)
library(gstat)
library(countrycode)
library(tidyr)
library(terra)

#load necessary data & build spatial objects
contdird <- read.csv("Simulations/contdird.csv") |> 
  filter(conttype == 1) |> 
  filter(year == 2016)

europe_boundaries <- vect("Simulations/Eurostat Shapefile")

codes <- c("AUS", "BEL", "BUL", "CZR", "DEN", "EST", "FIN", "FRN", "GMY", 
           "GRC", "HUN", "IRE", "ITA", "LAT", "LIT", "NTH", "NOR", "POL", "POR",
           "SLV", "SPN", "SWD", "SWZ", "UKG")

eur_boundaries_cropped <- crop(europe_boundaries, 
                               ext(2500000, 6500000, 1250000, 5416754.9853))
eur_boundaries_cropped$Country <- countrycode(unique(eur_boundaries_cropped$CNTR_ID),
                                              origin = "eurostat",
                                              destination = "cowc")
eur_boundaries_filtered <- subset(eur_boundaries_cropped, 
                                  eur_boundaries_cropped$Country %in% codes)

#Build neighborhood
nb <- build_neighborhood(eur_boundaries_filtered$Country, contdird)
