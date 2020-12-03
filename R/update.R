
#' @description Update a model to a new model by changing one formula
#' @param mod a HMM model object 
#' @param type are you updating the "hidden" or "obs" part of the model? 
#' @param i if "hidden" then i is the row of the structure containing the change
#'          if "obs" then i is the observation variable name
#' @param j if "hidden" then j is the column of the structure containing the change
#'          if "obs" then j is the parameter whose formula is to be changed 
#' @param change the change to make to the formula, see ?update.formula for details
#' @export
update.HMM <- function(mod, type, i, j, change) {
  
  if (type == "hidden") {
    # copy model components 
    new_obs <- mod$obs()$clone()
    copy_hid <- mod$hidden()$clone()
    # extract current structure 
    new_struct <- copy_hid$structure()
    # update relevant formula 
    new_struct[i, j] <- as.character(update(as.formula(new_struct[i, j]), change))
    # create new hidden sub-model component 
    new_hid <- MarkovChain$new(n_states = new_hid$nstates(), 
                               structure = new_struct, 
                               data = mod1$obs()$data()) 
  } else if (type == "obs") {
    # copy model components 
    copy_obs <- mod$obs()$clone()
    new_hid <- mod$hidden()$clone()
    # get formulas 
    forms <- copy_obs$formulas(raw = TRUE)
    # make change 
    forms[[i]][[j]] <- update(forms[[i]][[j]], change)
    # create new obs 
    new_obs <- Observation$new(data = copy_obs$data(), 
                               dists = copy_obs$dists(), 
                               n_states = copy_obs$nstates(), 
                               par = copy_obs$par(), 
                               formulas = forms)
  }
  # create new HMM object 
  new_mod <- HMM$new(obs = new_obs, hid = new_hid, init = mod)
  # fit new model 
  new_mod$fit() 
  # return fitted model 
  return(new_mod)
}
