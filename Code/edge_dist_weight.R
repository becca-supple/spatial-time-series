edge_dist_weight <- function(graphs, percent = TRUE, names = FALSE, weights = 1){
  
  n_rep <- length(weights)
  
  if(n_rep == 1){
    dist <- edge_dist(graphs = graphs, percent = percent, names = names)
  } else{
    
    
    n <- length(graphs) / n_rep
    
    if(ceiling(n) != floor(n)){
      stop("\nVector of weights must be a factor of total number of graphs.")
    }
    
    if(abs(sum(weights) - 1) > .Machine$double.eps^0.5){
      stop("\nWeights must sum to 1.")
    }
    
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

    # start_i <- as.numeric(str_split_1(names(graphs[1]), pattern = "_")[1])
    
    for(i in seq(1, n)){
      
      edges_i <- edge_dist(graphs[(i - 1)*n_rep + 1], percent = FALSE, names = FALSE)
      edges_i[,3:9] <- edges_i[,3:9] * weights[1]
      
      for(rep in seq(from = 2, to = n_rep)){
        
        edges_rep <- edge_dist(graphs[(i - 1)*n_rep + rep], percent = FALSE, names = FALSE)
        edges_rep[,3:9] <- edges_rep[,3:9] * weights[rep]
        
        edges_i[,3:9] <- edges_i[,3:9] + edges_rep[,3:9]
        
      }
      
      dist[,3:9] <- dist[,3:9] + edges_i[,3:9]
      
    }
    
    if(percent){
      dist[,3:9] <- (dist[,3:9]/n) * 100
    }
    
    if(names){
      
      #name the edge types
      edges <- list(E1 = "--",
                    E2 = "<-",
                    E3 = "<>",
                    E4 = "->",
                    E5 = "<*",
                    E6 = "*>")
      
      colnames(dist)[4:9] <- edges
    }

  }
  
  class(dist) <- c("edgedist", class(dist))
  return(dist)
  
}
