contains_truth <- function(graph_list,
                           true_edges = true_graph, 
                           allow_undirected = TRUE,
                           return_tprs = FALSE){
  
  same_as_truth <- rep(FALSE, length(graph_list))
  
  if(return_tprs){
    tprs <- numeric(length(graph_list))
  }
    
  force(true_edges)
  force(allow_undirected)
  
  for(g in seq(1, length(graph_list))){
    
    e <- edge_dist(graph_list[[g]])
    tpr <- accuracy_graph(e, true_edges = true_edges, allow_undirected = allow_undirected)
    
    if(tpr == 100){
      same_as_truth[g] <- TRUE
    }
    
    if(return_tprs){
      tprs[g] <- tpr
    }
    
  }
  
  p <- sum(same_as_truth)/length(same_as_truth)
  
  if(return_tprs){
    return(list(percent_containing_truth = p, tprs = tprs))
  }else{
    return(p)
  }
}
