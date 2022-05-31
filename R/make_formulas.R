
#' Process formulas and store in nested list
#' 
#' @param input_forms Nested list of formulas, with two levels: observed variable, 
#' and parameter of the observation distribution. The formulas can contain 
#' state-specific terms, e.g. "~ state1(x1) + x2".
#' @param var_names character vector name of each observation variable 
#' @param par_names list with element for each observation variable that
#' contains character vector of name of each parameter in its distribution 
#' @param n_states Number of states
#' 
#' @details Formulas for the observation parameters can be different for the
#' different states, using special functions of the form "state1", "state2", etc.
#' This method processes the list of formulas passed by the user to extract the 
#' state-specific formulas. Missing formulas are assumed to be intercept-only ~1. 
#' 
#' @return Nested list of formulas, with three levels: observed variable,
#' parameter of the observation distribution, and state.
#' 
#' @examples
#' input_forms <- list(step = list(shape = ~ state1(x1) + x2,
#'                                 scale = ~ x1),
#'                     count = list(lambda = ~ state1(x1) + state2(s(x2, bs = "cs"))))
#'
#' make_formulas(input_forms = input_forms, 
#'               var_names = names(input_forms), 
#'               par_names = lapply(input_forms, names), 
#'               n_states = 2)
make_formulas <- function(input_forms, 
                          var_names, 
                          par_names,
                          n_states) {
  # Output list
  output_forms <- list()
  
  # Get formula names
  form_names <- names(input_forms)
  
  # Loop over observed variables
  for(i in 1:length(var_names)) {
    # List of formulas for this variable, if any are given 
    mch <- match(var_names[i], form_names)
    if(!is.na(mch)) {
      var_forms <- input_forms[[mch]]
    } else {
      var_forms <- NULL
    }
    
    # Updated list of formulas, with extra level for state
    var_forms_new <- list()
    
    # Loop over parameters
    for(j in 1:length(par_names[[i]])) {
      # Formula for this parameter, if it is given 
      par_mch <- match(par_names[[i]][[j]], names(var_forms))
      if (!is.na(par_mch)) {
        form <- var_forms[[par_mch]]
      } else {
        form <- as.formula("~ 1")
      }
      
      # Terms object for this formula
      form_terms <- terms(form, specials = paste0("state", 1:n_states))
      
      # Term labels for this formula
      labs <- attr(form_terms, "term.labels")
      
      # Extract covariate names (remove special functions, e.g. state1, state2...)
      # The regular expression means the following:
      # ^: start of string
      # state[0-9]*: "state" followed by any number of characters between 0 and 9
      # \\(: opening bracket
      # (.*): grab what is between the brackets (then accessible as "\\1")
      # \\): closing bracket
      # $: end of string
      covs <- gsub(pattern = "^state[0-9]*\\((.*)\\)$", 
                   replacement = "\\1",  x = labs)
      
      # Find covariates which don't appear in any special function
      # (i.e. covariates that are included in all states)
      which_all_states <- which(!seq_along(labs) %in% unlist(attr(form_terms, "specials")))
      
      # Initialise list of state-specific formulas
      state_forms <- list()
      
      # Loop over states
      for(s in 1:n_states) {
        # Find covariates included in this state
        which_this_state <- attr(form_terms, "specials")[[paste0("state", s)]]
        
        # Initialise new formula
        new_form <- "~ 1"
        
        # Loop over terms that need to be added to the formula
        for(k in which_this_state)
          new_form <- paste0(new_form, " + ", covs[k])
        for(k in which_all_states)
          new_form <- paste0(new_form, " + ", covs[k])
        
        state_forms[[paste0("state", s)]] <- as.formula(new_form)
      }
      
      # Updated list of formulas for this parameter
      var_forms_new[[j]] <- state_forms
      names(var_forms_new)[j] <- par_names[[i]][[j]]
    }
    
    # Updated list of formulas for this variable
    output_forms[[i]] <- var_forms_new
    names(output_forms)[i] <- var_names[i]
  }
  
  return(output_forms)
}
