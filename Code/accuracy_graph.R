accuracy_graph <- function(est_edges, true_edges = true_graph, allow_undirected = "both"){
  
  # #name the edge types
  # edges <- list(E1 = "--",
  #               E2 = "<-",
  #               E3 = "<>",
  #               E4 = "->",
  #               E5 = "<*",
  #               E6 = "*>")
  
  if("edgedist" %in% class(est_edges) == FALSE){
    stop("est_edges must be an edgedist object")
  }
  
  n_edges <- nrow(est_edges)
  tps <- 0
  tps_d <- 0
  force(true_edges)
  
  for(e in seq(1, n_edges)){
    
    edge_t <- which.max(true_edges[e, 3:9])[[1]] + 2
    
    if(allow_undirected == TRUE || allow_undirected == "both"){
      
      if(edge_t == 8){
        cols <- c(4, 5, 6, 8)
      }else{
        if(edge_t == 9){
          cols <- c(4, 6, 7, 9)
        }else{
          cols <- edge_t
        }
      }
      
      edge_e <- sum(est_edges[e, cols])
      
      tps <- tps + edge_e
      
    }
    
    if(allow_undirected == FALSE || allow_undirected == "both"){
      
      edge_e <- est_edges[e, edge_t]
      tps_d <- tps_d + edge_e
      
    }
    
  }
  
  if(allow_undirected == TRUE){
    
    tpr <- tps/n_edges
    return(tpr)
    
  }else{
    
    if(allow_undirected == FALSE){
      
      tpr_d <- tps_d/n_edges
      return(tpr_d)
      
    }else{
      
      tpr <- tps/n_edges
      tpr_d <- tps_d/n_edges
      return(list(exact = tpr_d, undirected = tpr))
      
    }
  }
}
