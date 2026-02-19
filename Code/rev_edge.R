rev_edge <- function(edge){
  
  edge <- case_when(edge == "->" ~ "<-",
                    edge == "<-" ~ "->",
                    edge == "<*" ~ "*>",
                    edge == "*>" ~ "<*",
                    TRUE ~ edge)
  
  return(edge)
}
