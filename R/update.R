
#' Update a model to a new model by changing one formula
#' 
#' @param object HMM model object 
#' @param type Character string for the part of the model that is updated
#' (either "hid" or "obs") 
#' @param i If type = "hid" then i is the row of the formula containing the 
#' change. If type = "obs" then i is the observation variable name.
#' @param j If type = "hid" then j is the column of the formula containing the 
#' change. If type = "obs" then j is the parameter whose formula is to be 
#' changed. 
#' @param change The change to make to the formula, see ?update.formula for 
#' details.
#' @param fit If FALSE then change is made but model is not re-fit.
#' @param silent If TRUE then no model fitting output is given
#' @param ... Additional arguments are ignored (for compatibility with generic
#' S3 method)
#' 
#' @export
update.HMM <- function(object, type, i, j, change, fit = TRUE, 
                       silent = FALSE, ...) {
  
  # make sure data has state column added back in 
  dat <- object$obs()$data() 
  if (!all(is.na(object$obs()$known_states(mat = FALSE)))) {
    dat$state <- object$obs()$known_states(mat = FALSE)
  }
  
  if (type == "hid") {
    # copy model components 
    new_obs <- object$obs()$clone()
    copy_hid <- object$hid()$clone()
    # extract current formula 
    new_formula <- copy_hid$formula()
    # update relevant formula 
    new_formula[i, j] <- 
      as_character_formula(update(as.formula(new_formula[i, j]), change))
    # create new hidden sub-model component 
    new_hid <- MarkovChain$new(n_states = object$hid()$nstates(), 
                               formula = new_formula, 
                               data = dat, 
                               stationary = object$hid()$stationary()) 
  } else if (type == "obs") {
    # copy model components 
    copy_obs <- object$obs()$clone()
    new_hid <- object$hid()$clone()
    # get formulas 
    forms <- copy_obs$formulas(raw = TRUE)
    # make change 
    forms[[i]][[j]] <- update(forms[[i]][[j]], change)
    # create new obs 
    new_obs <- Observation$new(data = dat, 
                               dists = copy_obs$dists(), 
                               n_states = copy_obs$nstates(), 
                               par = copy_obs$inipar(), 
                               formulas = forms)
  }
  # create new HMM object 
  new_mod <- HMM$new(obs = new_obs, hid = new_hid, init = object, 
                     fixpar = object$fixpar())
  # fit new model 
  if(fit) new_mod$fit(silent = silent) 
  # return fitted model 
  return(new_mod)
}
