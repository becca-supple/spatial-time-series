edge_dist <- function(graphs, percent = TRUE, names = FALSE){
  
  #name the edge types
  edges <- list(E1 = "--",
                E2 = "<-",
                E3 = "<>",
                E4 = "->",
                E5 = "<*",
                E6 = "*>")
  
  #format as a list of graphs even if only 1 entered
  if(is.data.frame(graphs)){
    graphs <- list(graphs)
  }
  
  ###
  # Find distribution of  edges
  ###
  
  if (is.null(graphs[[1]])) {
    first_real <- which(sapply(graphs, is.null) == FALSE)[1]
    if(is.na(first_real)){
      return(NULL)
    }
    Xjs <- graphs[[first_real]]$Xj
    Xis <- graphs[[first_real]]$Xi
  } else {
    Xjs <- graphs[[1]]$Xj
    Xis <- graphs[[1]]$Xi
  }
  
  #set up dataframe of results
  dist <- data.frame(Xj = Xjs, #Xj column
                     Xi = Xis, #Xi column
                     Null = 0, #Initially 0 of each edge type
                     E1 = 0,
                     E2 = 0,
                     E3 = 0,
                     E4 = 0,
                     E5 = 0,
                     E6 = 0)
  
  #for each pair, find edge types in each graph
  for(pair in seq(1, nrow(dist))){
    
    Xj <- dist$Xj[pair]
    Xi <- dist$Xi[pair]
    
    for(g in seq(1, length(graphs))){
     
      graph <- graphs[[g]]
      
      if(is.null(graph)){
        next
      }
      
      edge <- graph[graph$Xj == Xj & graph$Xi == Xi, 3][[1]]
      
      #if an edge exists between this pair
      if(is.na(edge) == FALSE){
        
        #and increment the relevant column
        dist[[names(which(edges == edge))]][pair] <- dist[[names(which(edges == edge))]][pair] + 1
        
      }else{
        
        #if no edge, increment null column
        dist$Null[pair] <- dist$Null[pair] + 1
        
      }
    }
  }
  
  if(percent){
    #report values as percent of total
    n_graphs <- length(graphs) - sum(sapply(graphs, is.null))
    dist[,3:9] <- (dist[,3:9]/n_graphs) * 100
    
  }
  
  if(names){
    colnames(dist)[4:9] <- edges
  }
  
  class(dist) <- c("edgedist", class(dist))
  return(dist)
  
}
