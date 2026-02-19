accuracy_params <- function(graph_list, 
                       parameter_df = params, 
                       true_edges = true_graph, 
                       null_edges = true_null_graph,
                       n_reps = 2){
  
  df_list <- list()
  null_is <- which(params$strength == "none")
  
  for(varpair in seq(1, length(combn(colnames(params), 2)[1,]))){
    
    var1 <- combn(colnames(params), 2)[1,varpair]
    var2 <- combn(colnames(params), 2)[2,varpair]
    
    vals1 <- unique(params[[var1]])
    vals2 <- unique(params[[var2]])
    
    vardf <- expand.grid(vals1, vals2)
    colnames(vardf) <- c(var1, var2)
    vardf$tpr <- 0
    vardf[[var1]] <- as.factor(vardf[[var1]])
    vardf[[var2]] <- as.factor(vardf[[var2]])
    
    for(combo in seq(1, nrow(vardf))){
      
      is <- which(params[[var1]] == vardf[[var1]][combo] & 
                    params[[var2]] == vardf[[var2]][combo])
      
      tprs <- 0
      
      for(i in is){
        
        #cat(paste0("\n", i, ": ", seq((n_reps*(i-1) + 1), (n_reps*(i-1) + n_reps))))
        graphs <- graph_list[seq((n_reps*(i-1) + 1), (n_reps*(i-1) + n_reps))]
        edges <- edge_dist(graphs)
        
        if(i %in% null_is){
          tprs <- tprs + accuracy_graph(edges, true_edges = null_edges, allow_undirected = TRUE)
        }else{
          tprs <- tprs + accuracy_graph(edges, true_edges = true_edges, allow_undirected = TRUE)
        }
        
      }
      
      vardf$tpr[combo] <- sum(tprs)/(length(is))
      
    }
    
    df_list[[paste(var1, var2, sep = "-")]] <- vardf
    
  }
  
  class(df_list) <- c("tpr_df_list", class(df_list))
  return(df_list)
  
}
