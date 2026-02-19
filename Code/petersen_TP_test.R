# see Petersen, A. H. (2024). Are you doing better than random guessing? A call for using negative controls when evaluating causal discovery algorithms. Pre-print.

petersen_TP_test <- function(est_edges = NULL, true_edges = NULL, alpha = 0.05){
  
  if(is.null(est_edges)){
    error("Must provide (an) estimated graph(s) or edge distribution\n")
  }else{
    if("icgraph" %in% class(est_edges)){
    est_edges <- edge_dist(est_edges)
    }else{
      if("edgedist" %in% class(est_edges) == FALSE){
        error("Incorrect class for estimated graph(s) or edge distribution\n")
      }
    }
  }
  
  if(is.null(true_edges)){
    error("Must provide a true graph or edge distribution\n")
  }else{
    if("icgraph" %in% class(true_edges)){
      true_edges <- edge_dist(true_edges)
    }else{
      if("edgedist" %in% class(true_edges) == FALSE){
        error("Incorrect class for true graph or edge distribution\n")
      }
    }
  }
  
  true_adj <- max.col(true_edges[,3:9]) != 1
  est_adj <- max.col(est_edges[,3:9]) != 1
  
  m_max <- nrow(est_edges)
  m_true <- sum(true_adj)
  m_est <- sum(est_adj)
  
  TP_obs <- sum(est_adj & true_adj)
  
  p <- ifelse(TP_obs == 0, 1,
              1 - phyper(TP_obs - 1, m_max, m_true, m_est))
  
  list2return <- list(TP_obs = TP_obs, p = p, alpha = alpha, 
                      m_max = m_max, m_true = m_true, m_est = m_est)
  class(list2return) <- c("petersen_TP_test", class(list2return))
  return(list2return)
  
}

summary.petersen_TP_test <- function(test_output){
  
  asterisk <- ifelse(test_output$p < test_output$alpha, "*", "")
  
  cat(paste("True positive hypothesis test (Petersen, 2024) for a graph with:\n",
            test_output$m_est, "estimated edges versus", test_output$m_true, 
            "true edges, out of", test_output$m_max, 
            "possible edges.\nTest statistic (no. true positives):", 
            test_output$TP_obs, "\nP(TP >= TP observed):", 
            test_output$p, asterisk, "\n"))
  
}
