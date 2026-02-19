conditional_independencies <- function(data, 
                                       max_lag = 4, 
                                       factor,
                                       time = "Year",
                                       family = "gaussian",
                                       modeltype = "separated",
                                       neighborhood = nb,
                                       alpha = 0.05,
                                       silent = TRUE,
                                       skeleton_names = FALSE,
                                       skeleton_only = FALSE,
                                       collider_rule = "conservative",
                                       independence_test = "GCM",
                                       k0 = 7,
                                       nsim = 499){
  
  # where data is a dataframe of the time series with a column for 
  # the spatial factor variable and a column for the time unit, max_lag
  # is the maximum time lag tested, factor is the name of the factor variable
  # column in data, and time is the name of the time unit column in data. 
  # family can optionally be set to a list in the same order of the variables
  # in the dataframe with other options of families taken by mgcv:bam. Serial
  # dependence is set to TRUE for all variables by default, but can be instead 
  # provided as a vector of T/F values corresponding to the order of variables
  # in the dataset. Alpha is set to 0.05 by default for hypothesis tests, but 
  # can be changed/corrected. Set skeleton_names = TRUE to replace numerical variable
  # indices with variable names in the skeleton. Set skeleton_only = TRUE to 
  # skip edge orienting steps and just return skeleton. Can set collider_rule to
  # "conservative", "majority", or "none."
  
  ####
  # Force variables needed for later functions
  ####
  
  modeltype <- modeltype
  k0 <- k0
  nsim <- nsim
  max_lag <- max_lag
  if(independence_test %in% c("GCM", "WGCM") == FALSE){
    stop("\nError specifying independence test. Only GCM or WGCM independence tests supported.")
  }
  
  ####
  # Create variable set S & conditional independencies dataframe
  ####
  
  #Check for columns with NAs
  original_col <- ncol(data)
  data <- data |> 
    select_if(~ all(!is.na(.)))
  if(original_col != ncol(data)){
    stop("One or more columns containing NAs have been detected.")
  }
  
  if("design" %in% class(data) == FALSE){
    
    #Set of variables indexed
    X <- seq(1, sum(colnames(data) %in% c(factor, time) == FALSE))
    n_X <- length(X)
    
    design <- shift_to_design(data, 
                              site_col = factor, 
                              time_col = time, 
                              max_lag = max_lag)
    
  }else{
    
    design <- data
    
    #Set of variables indexed
    n_X <- (ncol(design) - 1) / (1 + (max_lag * 2))
    X <- seq(1, n_X)
    
  }
  
  #Each variable assigned a distribution
  if(length(family) == 1){
    distributions <- rep(family, (2*max_lag + 1)*n_X)
  }else{
    if(length(family) == n_X){
      distributions <- unlist(lapply(family, rep, (2*max_lag + 1)))
    }else{
      stop("\nError. Family list length does not match number of variables.")
    }
  }
  
  #Residual storing environment
  res_store <- new.env()
  
  #All present versions of variables
  X_t <- 1 + (1 + 2*max_lag)*(X - 1)
  
  #List of lagged sets for each Xj_t
  Bj_ts <- list()
  
  ###
  # Lagged Adjacencies
  ###
  
  for(j in seq(1:length(X_t))){
    
    #Define column in design relevant to Xj_t and its distribution
    
    Xj_t <- X_t[j]
    Xj_t_dist <- distributions[Xj_t]
    
    #Initialize lagged conditioning set
    Bj_t <- numeric(max_lag * n_X)
    ind <- 1
    
    for(tau in seq(1, max_lag)){
      for(i in X){
        Bj_t[ind] <- 1 + tau + (1 + 2*max_lag)*(i - 1)
        ind <- ind + 1
      }
    }
    
    #Dataframe to store independence test results
    store_ind <- data.frame(matrix(nrow = length(Bj_t),
                                   ncol = 4))
    
    colnames(store_ind) <- c("Xi_past", "Xj_t", "Imin", "remove")
    
    #Define test pairs, initialize test stat at Inf
    store_ind$Xj_t <- Xj_t
    store_ind$Xi_past <- Bj_t
    store_ind$Imin <- Inf
    store_ind$remove <- 0
    
    ####
    # Find lagged set
    ####
    
    #Initialize conditioning set size at 0
    p <- 0
    
    #While there are variables left to test
    while(length(Bj_t) - 1 >= p){
      
      #For each variable in the set
      for(Xi_past in Bj_t){
        
        #are there non-Xi_past vars in Bj_t?
        if(length(setdiff(Bj_t, Xi_past)) >= p){
          Bj_t_i <- setdiff(Bj_t, Xi_past)
        }else{
          next
        }
        
        #relevant index of independence storing dataframe
        i <- which(store_ind$Xi_past == Xi_past)
        
        #relevant distribution
        Xi_past_dist <- distributions[Xi_past]
        
        #Find the next largest conditioning set
        if(p == 0){
          #if empty set set model_on = "null"
          S <- "null"
        }else{
          #or define as first p values of Bj_t
          S <- Bj_t_i[1:p]
        }
        
        #Test if Xi_past _||_ Xj_t | S
        
        if(Xi_past %in% S || Xj_t %in% S){
          stop("Self-model in lagged testing")
        }
        
        if(!silent){cat("\nFitting ", colnames(design)[Xi_past], " ~ ", S)}
        
        #Fit and extract residuals for a null model of Xi_past
        Xi_past_S <- find_or_fit(design = design,
                                    response = Xi_past,
                                    model_on = S,
                                    factor = factor,
                                    distribution = Xi_past_dist,
                                    Env = res_store,
                                    nb = neighborhood,
                                    modeltype = modeltype)
        
        if(!silent){cat("\nFitting ", colnames(design)[Xj_t], " ~ ", S)}
        
        #Fit and extract residuals for Xj_t | S
        Xj_t_S <- find_or_fit(design = design,
                              response = Xj_t,
                              model_on = S,
                              factor = factor,
                              distribution = Xj_t_dist,
                              Env = res_store,
                              nb = neighborhood,
                              modeltype = modeltype)
        
        if(independence_test == "GCM"){
          test_result <-  gcm.test(resid.XonZ = Xi_past_S, 
                                   resid.YonZ =  Xj_t_S)
        }else{
          if(is.character(S)){
            zed <- as.numeric(design[[factor]])
          }else{
            zed <- design[,S]
          }
          
          test_result <- WGCM_fix(resid.XonZ = Xi_past_S, 
                                  resid.YonZ =  Xj_t_S,
                                  Z = zed,
                                  k0 = k0,
                                  nsim = nsim)
        }

        
        #Replace test statistic if lower than previously calculated
        if(test_result$test.statistic < store_ind$Imin[i]){
          store_ind$Imin[i]  <- test_result$test.statistic
        }
        
        #Mark for removal if p-value < alpha
        if(test_result$p.value > alpha){
          store_ind$remove[i]  <- 1
        }
      }
      
      #Remove those variables that were marked from Bj_t
      Bj_t <- setdiff(Bj_t, store_ind$Xi_past[store_ind$remove == 1])
      
      #store_ind <- store_ind[store_ind$Xi_past %in% Bj_t, ]
      
      #Sort the remaining set by Imin from largest to smallest
      Bj_t <- store_ind  |> 
        filter(Xi_past %in% Bj_t) |> 
        arrange(desc(Imin)) |> 
        pull(Xi_past)
      
      #Increment the size of the conditioning set
      p <- p + 1
      
    }
    
    #Store the finalized lagged set
    Bj_ts[[j]] <- Bj_t
    
  }
  
  ### 
  # Contemporaneous adjacencies
  ###
  
  #List of valid past indices
  pasts <- numeric(max_lag * n_X)
  ind <- 1
  
  for(tau in seq(1, max_lag)){
    for(i in X){
      pasts[ind] <- 1 + tau + (1 + 2*max_lag)*(i - 1)
      ind <- ind + 1
    }
  }
  
  #Initialize the graph
  G <- data.frame(Xj = c(rep(X_t, each = length(pasts)),
                         combn(X_t, 2)[1,]),
                  Xi = c(rep(pasts, n_X),
                         combn(X_t, 2)[2,]),
                  edge = NA)
  
  #Add links for all contemporaneous variables
  G <- G |> 
    mutate(edge = ifelse(Xi %in% X_t, "--", NA))
  
  #Add links from lagged adjacency step
  for(j in X){
    
    #if there were any
    if(length(Bj_ts[[j]]) > 0){
      
      #for each adjacency
      for(Xi in Bj_ts[[j]]){
        
        #add a marked edge from it to the present var in G
        G[G$Xj == X_t[j] & G$Xi == Xi, 3] <- "<*"
      }
    }
  }
  
  #Initialize all contemporaneous adjacency sets as all other present variables
  Aj_ts <- list()
  
  for(j in X){
    Aj_ts[[j]] <- X_t[X[-j]]
  }
  
  #Set up dataframe to store independence test results for all linked pairs
  store_ind <- G |> 
    filter(is.na(edge) == FALSE) |> 
    mutate(Imin = Inf,
           S = NA) |> 
    select(-edge)
  
  ###
  # Contemporaneous adjacency set testing
  ###
  
  p <- 0
  
  #While there is any linked pair for which the contemporaneous adjacency of Xj 
  #is at least 1 greater than the separating set size
  while(any(sapply(Aj_ts, function(x) {
    if (length(x) == 0) return(0)  # empty set treated as zero count
    length(x[x %in% G[!is.na(G$edge), 1]])
  }) - 1 >= p)) {
    
    js <- which(sapply(Aj_ts, function(x) {
      if(length(x) == 0) return(0)
      length(x[x %in% G[!is.na(G$edge), 1]])
    }) - 1 >= p)
    
    pairs <- which(is.na(G$edge) == FALSE & G$Xj %in% X_t[js])
    
    #For each pair in that condition
    for(pair in pairs){
      
      Xj_t <- G[pair, 1]
      j <- which(X_t == G[pair, 1])
      Xi_tau <- G[pair, 2]
      i <- ((Xi_tau - 1) %/% (1 + 2 * max_lag)) + 1
      
      #are there non-Xi_tau vars in Aj_t?
      Aj_t_i <- setdiff(Aj_ts[[j]], Xi_tau)
      
      #index in independence storing df
      ind <- which(store_ind$Xj == Xj_t & store_ind$Xi == Xi_tau)
      
      #Find all separating sets of size p in order
      if(p > 0){
        S_sets <- combn(as.list(Aj_t_i), p)
      }else{
        S_sets <- "null"
      }
      
      tested <- 0
      
      
      #While there are still untested sets
      while(tested < length(tested)){
        for(S in S_sets){
          
          #Set Z as the set containing S, lagged adj of Xj - Xi_tau, and lagged adj of Xi_tau
          Bj_t <- setdiff(Bj_ts[[j]], Xi_tau)
          
          #if tau = 0
          if(Xi_tau %in% X_t){
            Bi_t <- Bj_ts[[i]]
            
            #or if greater than 0
          }else{
            lag <- Xi_tau - X_t[i]
            Bi_t <- Bj_ts[[i]]  + lag
            
          }
          
          if(S == "null"){
            Z <- unique(c(Bj_t, Bi_t))
            if(length(Z) == 0){
              Z <- "null"
            }
          }else{
            Z <- unique(c(S, Bj_t, Bi_t))
          }
          
          #Test Xi_tau _||_ Xj_t | Z
          
          if(Xi_tau %in% Z  || Xj_t %in% Z){
            stop("Self-model in comtemp. testing")
          }
          
          if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}
          
          #Fit and extract residuals for a null model of Xi_past
          Xi_tau_Z <- find_or_fit(design = design,
                                     response = Xi_tau,
                                     model_on = Z,
                                     factor = factor,
                                     distribution = distributions[Xi_tau],
                                     Env = res_store,
                                     nb = neighborhood,
                                     modeltype = modeltype)
          
          if(!silent){cat("\nFitting ", colnames(design)[Xj_t], " ~ ", Z)}
          
          #Fit and extract residuals for Xj_t | Z
          Xj_t_Z <- find_or_fit(design = design,
                                response = Xj_t,
                                model_on = Z,
                                factor = factor,
                                distribution = distributions[Xj_t],
                                Env = res_store,
                                nb = neighborhood,
                                modeltype = modeltype)
          
          if(independence_test == "GCM"){
            test_result <-  gcm.test(resid.XonZ = Xi_tau_Z, 
                                     resid.YonZ =  Xj_t_Z)
          }else{
            if(is.character(Z)){
              zed <- as.numeric(design[[factor]])
            }else{
              zed <- design[,Z]
            }
            test_result <- WGCM_fix(resid.XonZ = Xi_tau_Z, 
                                    resid.YonZ =  Xj_t_Z,
                                    Z = zed,
                                    k0 = k0,
                                    nsim = nsim)
          }
          
          
          #Replace test statistic if lower than previously calculated
          if(test_result$test.statistic < store_ind$Imin[ind]){
            store_ind$Imin[ind]  <- test_result$test.statistic
          }
          
          #If p-value > alpha, delete link in G
          if(test_result$p.value > alpha){
            G[pair, 3] <- NA
            
            #and store S as the separating set
            store_ind$S[ind] <- paste(S, collapse = ";")
          }
          
          #Mark this set as tested
          tested <- tested + 1
          
          #Break if set was already found
          if(is.na(store_ind$S[ind]) == FALSE){
            tested <- length(tested)
            break
          }
          
        }
        
      }
      
      
    }
    
    #Increment p
    p <- p + 1
    
    #Reset contemporaneous sets based on the remaining links in G
    for(j in X){
      result <- store_ind |>
        filter(is.na(S)) |>
        filter(Xj == X_t[j] | Xi == X_t[j]) |>
        filter(Xi %in% X_t & Xj %in% X_t) |>
        arrange(desc(Imin)) |>
        mutate(other = ifelse(Xj == X_t[j], Xi, Xj)) |>
        pull(other)
      
      #Remove self (response) explicitly
      result <- setdiff(result, X_t[j])
      
      if(length(result) == 0){
        Aj_ts[[j]] <- numeric(0)  # assign empty numeric vector if no result
      }else{
        Aj_ts[[j]] <- result
      }
    }
    
  }
  
  ###
  # Format and return if you just want skeleton
  ###

  if(skeleton_only == TRUE){
    if(skeleton_names == TRUE){
      #replace indices with names for interpretability
      G <- G |> 
        mutate(Xj = colnames(design)[Xj],
               Xi = colnames(design)[Xi])
    }
    
    list2return <- list(skeleton = G, independencies = store_ind)
    return(list2return)
  }
  
  G_old <- G #store skeleton for later
  
  ###
  # Otherwise, ID unshielded triples and test for colliders
  ###
  
  #Set up df to store indices of unsh triples
  unsh_triples <- data.frame(Xi_tau = NA, Xk_t = NA, Xj_t = NA, ambiguous = NA,
                             flag = 1)
  
  testing <- TRUE
  
  #While there are unshielded triples to test
  while(testing == TRUE){
    
    ###
    # ID unshielded triples
    ###
    
    contemp_adj <- G |> 
      filter(edge == "--") |> 
      filter(Xi %in% X_t & Xj %in% X_t) 
    
    if(nrow(contemp_adj) < 1){
      unsh_triples <- unsh_triples[-1,]
      testing <- FALSE
      break
    }
    
    for(row in seq(1, nrow(contemp_adj))){
      
      #trying with this orientation of Xj/Xi first, then will do opposite
      Xk_t <- contemp_adj$Xj[row]
      Xj_t <- contemp_adj$Xi[row]
      
      #Is there any Xi_tau that has an arrow or edge into Xk_t but not adj to Xj_t?
      
      Xi_tau <- G |> #arrow or edge into Xk_t?
        filter(Xj == Xk_t | Xi == Xk_t) |> 
        filter(Xj != Xj_t & Xi != Xj_t) |> 
        filter(edge %in% c("--", "<*")) |> #unoriented edge if present, <* if from past
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)
      
      if(length(Xi_tau) == 0){
        next
      }
      
      Xi_tau <- G |> #no arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |> 
        filter(is.na(edge)) |> #only those with no edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)
      
      if(length(Xi_tau) == 0){
        next
      } 
      
      #there are instance(s) of an unshielded triple, add to df
      for(Xi in Xi_tau){
        
        if(is.na(unsh_triples[1,1])){ #replace 1st row if need be
          ind <- 1
        }else{
          ind <- nrow(unsh_triples) + 1 #otherwise add to end of df
        }
        
        unsh_triples[ind,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0
        )
      }
      
      #Reattempt with alternative Xk/Xj assignment
      Xk_t <- contemp_adj$Xi[row]
      Xj_t <- contemp_adj$Xj[row]
      
      #Is there any Xi_tau that has an arrow or edge into Xk_t but not adj to Xj_t?
      
      Xi_tau <- G |> #arrow or edge into Xk_t?
        filter(Xj == Xk_t | Xi == Xk_t) |> 
        filter(Xj != Xj_t & Xi != Xj_t) |> 
        filter(edge %in% c("--", "<*")) |> #unoriented edge if present, <* if from past
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)
      
      if(length(Xi_tau) == 0){
        next
      }
      
      Xi_tau <- G |> #no arrow or edge into Xj_t?
        filter(Xj == Xk_t | Xi == Xk_t) |> 
        filter(is.na(edge)) |> #only those with no edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)
      
      if(length(Xi_tau) == 0){
        next
      } 
      
      #there are instance(s) of an unshielded triple, add to df
      for(Xi in Xi_tau){
        
        if(is.na(unsh_triples[1,1])){ #replace 1st row if need be
          ind <- 1
        }else{
          ind <- nrow(unsh_triples) + 1 #otherwise add to end of df
        }
        
        unsh_triples[ind,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0
        )
      }
    }
    
    #if no unsh_triples were found at all, break
    if(is.na(unsh_triples[1,1])){ 
      testing <- FALSE
      break
    }
    
    ###
    # Test unshielded triples for colliders
    ###
    
    for(row in seq(1, nrow(unsh_triples))){
      
      Xi_tau <- unsh_triples$Xi_tau[row]
      i <- ((Xi_tau - 1) %/% (1 + 2 * max_lag)) + 1
      
      Xk_t <- unsh_triples$Xk_t[row]
      
      Xj_t <- unsh_triples$Xj_t[row]
      j <- which(X_t == Xj_t)
      
      if(collider_rule == "none"){
        
        sep_set <- store_ind |> 
          filter((Xi == Xi_tau & Xj == Xj_t) | (Xi == Xj_t & Xj == Xi_tau)) |> 
          pull(S)
        
        if(length(sep_set) > 0){
          #If Xk_t is not in the separating set
          if(Xk_t %in% as.numeric(str_split(sep_set, ";", simplify = TRUE)) == FALSE){
            
            #Orient Xj_t -> Xk_t
            G <- G |> 
              mutate(
                edge = case_when(
                  Xj == Xj_t & Xi == Xk_t ~ "*>",
                  Xi == Xk_t & Xj == Xj_t ~ "<*",
                  TRUE ~ edge  # retain existing value
                )
              )
            
            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1
          }
        }
        
        if(row == nrow(unsh_triples)){
          break
        }
        next
      }
      
      Aj_t <- Aj_ts[[j]]
      if(Xi_tau %in% X_t){ #include contemp set of i if relevant
        Ai_t <- Aj_ts[[i]]
        C_set <- setdiff(unique(c(Aj_t, Ai_t)), Xi_tau)
      }else{
        C_set <- setdiff(Aj_t, Xi_tau)
      }
      
      C_set <- as.vector(C_set)
      
      #All subsets
      n <- length(C_set)
      S_s <- list()
      
      if (n > 0) {
        for (w in 1:(2^n - 1)) {  # skip 0 to avoid empty set
          mask <- as.logical(intToBits(w))[1:n]
          S_s[[length(S_s) + 1]] <- C_set[mask]
        }
      }
      
      for(S in S_s){
        
        #Set Z as the set containing S, lagged adj of Xj - Xi_tau, and lagged adj of Xi_tau
        Bj_t <- setdiff(Bj_ts[[j]], Xi_tau)
        
        #if tau = 0
        if(Xi_tau %in% X_t){
          Bi_t <- Bj_ts[[i]]
          #or if greater than 0
        }else{
          lag <- Xi_tau - X_t[i]
          Bi_t <- Bj_ts[[i]]  + lag
        }
        
        Z <- unique(c(S, Bj_t, Bi_t))
        
        if(Xi_tau %in% Z || Xj_t %in% Z){
          stop("Self model in collider testing")
        }
        
        if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}
        
        #Test Xi_tau _||_ Xj_tau | Z:
        #Fit and extract residuals for a model of Xi_past
        Xi_tau_Z <- find_or_fit(design = design,
                                   response = Xi_tau,
                                   model_on = Z,
                                   factor = factor,
                                   distribution = distributions[Xi_tau],
                                   Env = res_store,
                                   nb = neighborhood,
                                   modeltype = modeltype)
        
        if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}
        
        #Fit and extract residuals for Xj_t | Z
        Xj_t_Z <- find_or_fit(design = design,
                              response = Xj_t,
                              model_on = Z,
                              factor = factor,
                              distribution = distributions[Xj_t],
                              Env = res_store,
                              nb = neighborhood,
                              modeltype = modeltype)

        if(independence_test == "GCM"){
          test_result <-  gcm.test(resid.XonZ = Xi_tau_Z, 
                                   resid.YonZ =  Xj_t_Z)
        }else{
          if(is.character(Z)){
            zed <- as.numeric(design[[factor]])
          }else{
            zed <- design[,Z]
          }
          test_result <- WGCM_fix(resid.XonZ = Xi_tau_Z, 
                                  resid.YonZ =  Xj_t_Z,
                                  Z = zed,
                                  k0 = k0,
                                  nsim = nsim)
        }
        
        if(test_result$p.value > alpha){ #if independent
          
          #store as sep. set
          store_ind[nrow(store_ind) + 1,] <- c(
            #Xj
            Xj_t,
            #Xi
            Xi_tau,
            #Imin
            test_result$test.statistic,
            #S
            paste(S, collapse = ";")
          )
          
        }
        
      }
      
      #list of all separating sets
      all_sets <- store_ind |> 
        filter(Xj == Xj_t & Xi == Xi_tau) |> 
        pull(S)
      
      #If no separating sets found, mark as ambiguous:
      if(length(all_sets) == 0){
        unsh_triples$ambiguous[row] <- 1
      }else{ #If some found:
        
        n_k <- mean(grepl(paste0("(?<!\\d)", Xk_t, "(?!\\d)"), 
                          unique(all_sets), perl=TRUE))
        
        if(collider_rule == "conservative"){
          
          if(n_k == 0){
            #Orient Xj_t -> Xk_t
            G <- G |> 
              mutate(
                edge = case_when(
                  Xj == Xj_t & Xi == Xk_t ~ "*>",
                  Xi == Xk_t & Xj == Xj_t ~ "<*",
                  TRUE ~ edge  # retain existing value
                )
              )
            
            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1
            
          }else{
            if(n_k == 1){
              #leave as is
              next
            }else{
              #if any other value, mark as ambiguous
              unsh_triples$ambiguous[row] <- 1
            }
          }
        }else{ #if collider_rule is majority or unset
          if(n_k < 0.5){
            #Orient Xj_t -> Xk_t
            G <- G |> 
              mutate(
                edge = case_when(
                  Xj == Xj_t & Xi == Xk_t ~ "*>",
                  Xi == Xk_t & Xj == Xj_t ~ "<*",
                  TRUE ~ edge  # retain existing value
                )
              )
            
            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1
            
          }else{
            if(n_k > 0.5){
              #leave as is
              next
            }else{
              #if = 0.5 or NA, mark as ambiguous
              unsh_triples$ambiguous[row] <- 1
            }
          }
        }
      }
    }
    
    #Once all unsh_triples tested, remove those flagged
    unsh_triples <- unsh_triples |> 
      filter(flag == 0)
    
    testing <- FALSE
  }
  
  ###
  # Orienting unambigious triples
  ###
  
  if(skeleton_names == TRUE){
    #replace indices with names for interpretability
    G <- G |> 
      mutate(Xj = colnames(design)[Xj],
             Xi = colnames(design)[Xi])
    
    list2return <- list(skeleton = G, independencies = store_ind, triples = unsh_triples)
  }else{
    list2return <- list(skeleton = G, independencies = store_ind, triples = unsh_triples, names = colnames(design))
  }
  
  return(list2return)
  
}
