build_neighborhood <- function(site_names, #names of sites you're using
                               contiguity_df, #data frame with all site combinations and a column that == 1 if contiguous
                               site1 = "state1ab", #name of column with site 1
                               site2 = "state2ab", #name of column with site 2
                               contiguity_col = "conttype" #name of column with contiguity info
){
  
  nb <- list()
  
  #Filter to just contiguities
  contiguity_df <- contiguity_df[contiguity_df[[contiguity_col]] == 1,]
  
  for(c in site_names){
    
    #find all contiguous other sites
    connections <- contiguity_df |> 
      mutate(connections = ifelse(contiguity_df[[site1]] == c, paste(contiguity_df[[site2]]),
                                  ifelse(contiguity_df[[site2]] == c, paste(contiguity_df[[site1]]), NA))) |> 
      select(connections)
    
    #filter to sites we're interested in
    connections <- unique(connections$connections[is.na(connections$connections) == F])
    connections <- connections[connections %in% site_names]
    
    nb[[c]] <- connections
  }
  
  return(nb)
  
}
