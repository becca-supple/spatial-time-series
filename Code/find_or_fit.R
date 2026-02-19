find_or_fit <- function(design, #design matrix
                        response = a, #variable to fit
                        model_on = b, #variable(s) to predict with or "null" for null model
                        factor = factor, #name of factor variable
                        distribution = distributions[a], #distribution of a
                        Env = res_store,
                        nb = nb,
                        modeltype = modeltype
){
  
  cond_key <- paste(model_on, collapse = ";")
  key <- paste(response, cond_key, sep = "_")
  
  if(exists(key, envir = Env)){
    
    #extract from storage
    return(Env[[key]])
    
  }else{ #otherwise fit the and store residuals for later
    
    #force all arguments to pass to fit_custom_gam()
    design <- design
    model_on <- model_on
    factor <- factor
    distribution <- force(distribution)
    nb <- force(nb)
    modeltype <- modeltype
    
    #fit and extract residuals
    Env[[key]] <- fit_custom_gam(design, response, model_on, factor, 
                                 distribution, nb, modeltype)
    
    return(Env[[key]])
    
  }
  
}
