shift_to_design <- function(data, site_col = "Country", time_col = "Year", max_lag = 2){
  
  #Create a version of data with shifted versions of all time series
  n <- sum(data[[site_col]] == unique(data[[site_col]][1]))
  num_rows <- n - 2 * max_lag #number of rows in design per site_col level
  n_X <- ncol(data) - 2
  X <- seq(1, n_X)
  
  #make sure time and site columns are at end
  data <- data |> 
    select(-c(all_of(site_col), all_of(time_col)), c(all_of(site_col), all_of(time_col)))
  
  design <- data.frame(matrix(
    nrow = (num_rows)*length(unique(data[[site_col]])), 
    ncol = (2*max_lag + 1)*n_X + 1))
  
  colnames(design)[(2*max_lag + 1)*n_X + 1] <- site_col
  col_index <- 0
  
  #Name the columns
  for(d in X){
    
    variable <- colnames(data)[d]
    for(k in seq(0, 2*max_lag)){
      if(k > 0){
        colnames(design)[col_index + k + 1] <- paste(variable, "_tminus", k, sep = "")
      } else{
        colnames(design)[col_index + k + 1] <- variable
      }
    }
    col_index <- col_index + (2*max_lag + 1)
  }
  
  #Fill in data values
  start <- 1 
  for(fac in unique(data[[site_col]])){
    filtered <- dplyr::filter(data, .data[[site_col]] == fac) #filter time series
    col_index <- 0
    
    for(d in X){
      for(k in seq(0, 2*max_lag)){ #fill in data
        design[start:(start + num_rows - 1), col_index + k + 1] <- filtered[seq(
          2*max_lag + 1 - k, n - k), d]
      }
      col_index <- col_index + (2*max_lag + 1)
    }
    
    #label by factor value
    design[[site_col]][start:(start + num_rows - 1)] <- rep(fac, num_rows)
    
    #shift the data entry to account for previous entered values
    start <- start + num_rows 
  }
  
  #site_col variable as a factor
  design[[site_col]] <- factor(design[[site_col]])
  
  #Set class 
  class(design) <- c("design", class(design))
  
  return(design)
  
}
