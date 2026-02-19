fit_custom_gam <- function(design = design, #design matrix
                           response = a, #index in design of variable to fit
                           model_on = b, #indices in design of variable(s) to predict with or "null" for null model
                           factor = factor, #name of factor variable
                           distribution, #distribution of a
                           nb, #neighborhood structure
                           modeltype #separated vs tensor vs both
){
  
  nb <- force(nb)
  if (response %in% model_on) {
    stop("Self-model detected: response appears in predictors")
  }
  
  if("null" %in% model_on){
    
    #Create MRF smooth
    formula <- paste0("s(", factor, ", bs = 'mrf', xt = list(nb = nb))") 
    #Connect to response
    formula <- paste(colnames(design)[response], formula, sep = " ~ ")
    #Format
    formula <- as.formula(formula)
    
    fam_fun <- get(distribution)
    
    null_model <- suppressWarnings(gam(formula, family = fam_fun(), data = design, 
                      method = "REML"))
    res <- residuals(null_model, type = "response")
    
    return(res)
    
  }else{
    
    if(modeltype == "separated"){
      #Initialize formula
      formula <- paste0(colnames(design)[response], " ~ ")
      
      
      #Find the right column in design for each predictor
      vars_c_design <- colnames(design)[model_on]
      
      for(c_ind in seq(1, length(model_on))){
        
        #Combine into a single formula
        
        c_formula <- paste("s(", vars_c_design[c_ind], ", bs = 'tp')", sep = "")
        
        if(c_ind == 1){
          formula <- paste(formula, c_formula)
        }else{
          formula <- paste(formula, c_formula, sep = " + ")
        }
        
      }
      
      #Add spatial smooth
      formula <- paste(formula, 
                       paste0("s(", factor, ", bs = 'mrf', xt = list(nb = nb))"),
                       sep = " + ")
      
      formula <- as.formula(formula)
      
      fam_fun <- get(distribution)
      
      #Fit GAM
      fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(), 
                                      data = design, method = "REML"))
      
      #Return residuals
      res <- residuals(fit_bam, type = "response")
      
      return(res)
    }else{
      if(modeltype == "tensor"){
        
        #Initialize formula
        formula <- paste0(colnames(design)[response], " ~ ")
        
        
        #Find the right column in design for each predictor
        vars_c_design <- colnames(design)[model_on]
        
        for(c_ind in seq(1, length(model_on))){
          
          #Combine into a single formula
          
          c_formula <- paste0("te(", vars_c_design[c_ind], 
                             ",", 
                             factor, 
                             ", bs = c('tp', 'mrf'), xt = list(nb = nb))")
          
          if(c_ind == 1){
            formula <- paste(formula, c_formula)
          }else{
            formula <- paste(formula, c_formula, sep = " + ")
          }
          
        }
        
        formula <- as.formula(formula)
        
        fam_fun <- get(distribution)
        
        #Fit GAM
        fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(), 
                                        data = design, method = "REML"))
        
        #Return residuals
        res <- residuals(fit_bam, type = "response")
        
      }else{
        if(modeltype == "both"){
          
          #Initialize formula
          formula <- paste0(colnames(design)[response], " ~ ")
          
          #Find the right column in design for each predictor
          vars_c_design <- colnames(design)[model_on]
          
          for(c_ind in seq(1, length(model_on))){
            
            #Combine into a single formula
            
            c_formula <- paste0("s(", vars_c_design[c_ind], ", bs = 'tp') + ", 
                                "te(", vars_c_design[c_ind], 
                                ",", 
                                factor, 
                                ", bs = c('tp', 'mrf'), xt = list(nb = nb))")
            
            if(c_ind == 1){
              formula <- paste(formula, c_formula)
            }else{
              formula <- paste(formula, c_formula, sep = " + ")
            }
            
          }
          
          #Add spatial smooth
          formula <- paste(formula, 
                           paste0("s(", factor, ", bs = 'mrf', xt = list(nb = nb))"),
                           sep = " + ")
          
          formula <- as.formula(formula)
          
          fam_fun <- get(distribution)
          
          #Fit GAM
          fit_bam <- suppressWarnings(gam(formula = formula, family = fam_fun(), 
                                          data = design, method = "REML"))

          
          #Return residuals
          res <- residuals(fit_bam, type = "response")
          
        }else{
          stop(paste0("\nUnknown model type ", modeltype,
                      ". Options are 'separated', 'tensor', or 'both'."))}
      }
    }
  }
}
