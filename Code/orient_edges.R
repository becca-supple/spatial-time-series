orient_edges <- function(independencies, #outputs from conditional independence algorithm
                         max_lag = 2, #max_lag for index calc purposes
                         names = TRUE #replace indices with variable names
){
  
  G <- independencies$skeleton
  trips <- independencies$triples
  X <- seq(1, (length(independencies$names) - 1)/ (2 * max_lag + 1))
  X_t <- 1 + (1 + 2*max_lag)*(X - 1)
  
  ###
  # ID unambiguous triples
  ###
  
  #Start with unambiguous triples ID'd in last algorithm
  
  unamb_trips <- trips |> 
    filter(ambiguous == 0) |> 
    mutate(cycle = 0)
  
  #Add triples that form a potential cycle
  
  contemp_adj <- G |> 
    filter(is.na(edge) == FALSE) |> 
    filter(Xi %in% X_t & Xj %in% X_t)
  
  if(nrow(contemp_adj) > 0){ #if there are contemporaneous adjacencies
    
    for(row in seq(1, nrow(contemp_adj))){
      
      #trying with this orientation of Xj/Xi first, then will do opposite
      Xk_t <- contemp_adj$Xj[row]
      Xj_t <- contemp_adj$Xi[row]
      
      #Is there any Xi_tau that has an arrow or edge into Xk_t and adj to Xj_t?
      
      Xi_tau <- G |> #arrow or edge into Xk_t?
        filter(Xj == Xk_t | Xi == Xk_t) |> 
        filter(is.na(edge) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)
      
      Xi_tau <- G |> #arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |> 
        filter(is.na(edge) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)
      
      if(length(Xi_tau) == 0){
        next
      }
      
      #if there are relevant Xi_tau, add to the triples df
      for(Xi in Xi_tau){
        
        unamb_trips[nrow(unamb_trips) + 1,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0,
          #cycle
          1
        )
      }
    }
  }
  
  while(nrow(unamb_trips) > 0){
    
    for(tr in seq(1, nrow(unamb_trips))){
      
      Xi_tau <- unamb_trips$Xi_tau[tr]
      Xk_t <- unamb_trips$Xk_t[tr]
      Xj_t <- unamb_trips$Xj_t[tr]

      edgeXi_tau_Xk_t <- G |> 
        filter(Xj == Xi_tau & Xi == Xk_t |
                 Xj == Xk_t & Xi == Xi_tau) |> 
        mutate(print_edge = ifelse(Xj == Xi_tau, edge, rev_edge(edge))) |> 
        pull(print_edge)
      
      edgeXk_t_Xj_t <- G |> 
        filter(Xj == Xk_t & Xi == Xj_t |
                 Xj == Xj_t & Xi == Xk_t) |> 
        mutate(print_edge = ifelse(Xj == Xk_t, edge, rev_edge(edge))) |> 
        pull(print_edge)
      
      ###
      # Rule 1: Non-colliders with an arrow into the middle must be chains
      ###
      
      if(edgeXi_tau_Xk_t %in% c("->", "*>") && 
         edgeXk_t_Xj_t == "--" && 
         unamb_trips$cycle[tr] == 0){
        
        #If an arrow into Xk_t from Xi_tau; orient from Xk_t to Xj_t
        
        G <- G |> 
          mutate(edge = case_when(
            Xj == Xk_t & Xi == Xj_t ~ "*>",
            Xj == Xj_t & Xi == Xk_t ~ "<*",
            TRUE ~ edge
          ))
        
        #Flag for removal 
        
        unamb_trips$flag[tr] <- 1
        
      }
      
      ###
      # Rule 2: Avoid cycles
      ###
      
      if(edgeXi_tau_Xk_t == "*>" && 
         edgeXk_t_Xj_t == "*>" && 
         unamb_trips$cycle[tr] == 1){
        

        #If there is a chain with no arrow into Xj_t, orient Xi_tau -> Xj_t
        edgeXi_tau_Xj_t <- G |> 
          filter(Xj == Xi_tau & Xi == Xj_t |
                   Xj == Xj_t & Xi == Xi_tau) |> 
          mutate(print_edge = ifelse(Xj == Xi_tau, edge, rev_edge(edge))) |> 
          pull(print_edge)
        
        if(edgeXi_tau_Xj_t %in% c("--", "<-")){
          
          if(edgeXi_tau_Xj_t == "--"){
            G <- G |> 
              mutate(edge = case_when(
                Xj == Xi_tau & Xi == Xj_t & Xi_tau %in% X_t ~ "->",
                Xj == Xi_tau & Xi == Xj_t & Xi_tau %in% X_t == FALSE ~ "*>", #mark if from past
                Xj == Xj_t & Xi == Xi_tau & Xi_tau %in% X_t ~ "<-",
                Xj == Xj_t & Xi == Xi_tau & Xi_tau %in% X_t == FALSE ~ "<*",
                TRUE ~ edge
              ))
          }else{
            G <- G |> 
              mutate(edge = case_when(
                Xj == Xi_tau & Xi == Xj_t & Xi_tau %in% X_t ~ "<>",
                Xj == Xi_tau & Xi == Xj_t & Xi_tau %in% X_t == FALSE ~ "*>", #mark if from past
                Xj == Xj_t & Xi == Xi_tau & Xi_tau %in% X_t ~ "<>",
                Xj == Xj_t & Xi == Xi_tau & Xi_tau %in% X_t == FALSE ~ "<*",
                TRUE ~ edge
              ))
          }
        }
        
        #Flag for removal 
        
        unamb_trips$flag[tr] <- 1
        
      }
      
      ###
      # Rule 3: Orient central edges in double cycles
      ###
      
      if(edgeXi_tau_Xk_t == "--" && 
         edgeXk_t_Xj_t %in% c("->", "*>") && 
         unamb_trips$cycle[tr] == 1){
        
        #are there adjacent cycles we can use?
        Xl_t <- unamb_trips |> 
          filter(Xi_tau %in% c(Xi_tau, Xj_t) &
                   Xj_t %in% c(Xi_tau, Xj_t) &
                   Xk_t != Xk_t) |> 
          pull(Xk_t)
        
        Xl_t <- intersect(Xl_t, X_t) #Xl_t in the present
        
        if(length(Xl_t) > 0){ #if there are relevant Xls, do they have the right edges
          
          for(Xl in Xl_t){
            
            edgeXi_tau_Xl_t <- G |> 
              filter(Xj == Xi_tau & Xi == Xl |
                       Xj == Xl & Xi == Xi_tau) |> 
              mutate(print_edge = ifelse(Xj == Xi_tau, edge, rev_edge(edge))) |> 
              pull(print_edge)
            
            edgeXl_t_Xj_t <- G |> 
              filter(Xj == Xl & Xi == Xj_t |
                       Xj == Xj_t & Xi == Xl) |> 
              mutate(print_edge = ifelse(Xj == Xl, edge, rev_edge(edge))) |> 
              pull(print_edge)
            
            if(edgeXi_tau_Xl_t == "--" && 
               edgeXl_t_Xj_t %in% c("->", "*>")){ #If the triple has the right structure
              
              #Orient Xi_tau *> Xj_t
              G <- G |> 
                mutate(edge = case_when(
                  Xj == Xi_tau & Xi == Xj_t ~ "*>",
                  Xj == Xj_t & Xi == Xi_tau ~ "<*",
                  TRUE ~ edge
                ))
              
              #Flag for removal 
              
              unamb_trips$flag[tr] <- 1
            }
          }
        }
      }
    }

    #Remove any flagged entries
    unamb_trips_new <- unamb_trips |> 
      filter(flag == 0)
    
    #If no change in the df, break the loop
    if(identical(unamb_trips_new, unamb_trips)){
      unamb_trips <- matrix(nrow = 0, ncol = 1)
    }else{ #if changes, update to newer version
      unamb_trips <- unamb_trips_new
    }
  }
  
  ###
  # Format and return
  ###
  
  #Replace indices with variable names
  if(names){
    G <- G |> 
      mutate(Xj = independencies$names[Xj],
             Xi = independencies$names[Xi])
  }
  
  class(G) <- c("icgraph", class(G))
  
  return(G)
  
}
