
#' Process formulas and store in nested list
#' 
#' @param input_forms Nested list of formulas, with two levels: observed variable, 
#' and parameter of the observation distribution. The formulas can contain 
#' state-specific terms, e.g. "~ state1(x1) + x2".
#' @param n_states Number of states
#' 
#' @return Nested list of formulas, with three levels: observed variable,
#' parameter of the observation distribution, and state.
#' 
#' @example 
#' input_forms <- list(step = list(shape = ~ state1(x1) + x2,
#'                                 scale = ~ x1),
#'                     count = list(lambda = ~ state1(x1) + state2(s(x2, bs = "cs"))))
#'
#' make_formulas(input_forms = input_forms, n_states = 2)
make_formulas <- function(input_forms, n_states) {
  # Output list
  output_forms <- list()
  
  # Loop over observed variables
  for(i in 1:length(input_forms)) {
    # List of formulas for this variable
    var_forms <- input_forms[[i]]
    
    # Updated list of formulas, with extra level for state
    var_forms_new <- list()
    
    # Loop over parameters
    for(j in 1:length(var_forms)) {
      # Formula for this parameter
      form <- var_forms[[j]]
      
      # Terms object for this formula
      form_terms <- terms(form, specials = paste0("state", 1:n_states))
      
      # Term labels for this formula
      labs <- attr(form_terms, "term.labels")
      n_terms <- length(labs)
      
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
      which_all_states <- which(! (1:n_terms) %in% unlist(attr(form_terms, "specials")))
      
      # Loop over states
      state_forms <- list()
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
      names(var_forms_new)[j] <- names(var_forms)[j]
    }
    
    # Updated list of formulas for this variable
    output_forms[[i]] <- var_forms_new
    names(output_forms)[i] <- names(input_forms)[i]
  }
  
  return(output_forms)
}
