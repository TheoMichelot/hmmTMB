
#' Create model matrices
#' 
#' @param formulas List of formulas (possibly nested, e.g. for use within Observation)
#' @param data Data frame including covariates
#' @param new_data Optional new data set, including covariates for which
#' the design matrices should be created. This needs to be passed in addition
#' to the argument '\code{data}', for cases where smooth terms or factor
#' covariates are included, and the original data set is needed to determine
#' the full range of covariate values.
#' @param gam_args Named list of additional arguments for \code{mgcv::gam()},
#' such as knots.
#' 
#' @return A list of
#' \itemize{
#'   \item X_fe Design matrix for fixed effects
#'   \item X_re Design matrix for random effects
#'   \item S Smoothness matrix
#'   \item log_det_S Vector of log-determinants of smoothness matrices, one
#'   for each penalty
#'   \item ncol_fe Number of columns of X_fe for each parameter
#'   \item ncol_re Number of columns of X_re and S for each random effect
#'   \item L Matrix with log(lambda) = L * theta, mapping the smoothing
#'   parameters of each smooth to the weights of its penalties (mgcv's
#'   convention). It is the identity unless a smooth combines several
#'   penalties through fewer parameters.
#'   \item sp_names Name of each smoothing parameter
#'   \item theta_start Starting value of each log smoothing parameter
#' }
#'
#' @details
#' Note the two indices: \code{ncol_re} and \code{log_det_S} have one entry per
#' \emph{penalty}, while \code{L}, \code{sp_names} and \code{theta_start} have
#' one per \emph{smoothing parameter}. The two coincide for an ordinary smooth,
#' which has one of each.
#' 
#' @importFrom stats update predict
make_matrices = function(formulas, data, new_data = NULL, gam_args = NULL) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  names_fe <- NULL
  names_re <- NULL
  names_ncol_re <- NULL
  log_det_S <- NULL
  L_list <- list()
  sp_names <- NULL
  theta_start <- NULL
  start <- 1
  
  # Unlist formulas so that this function works both for Observation and MarkovChain
  forms <- unlist(formulas)
  names <- names(forms)
  
  # Loop over formulas
  for(k in seq_along(forms)) {
    form <- forms[[k]]
    
    # Check that tensor products aren't used (not supported)
    check_smooths(form)
    
    # Check that random effect variables are factors
    var_re <- find_re(form)
    for(var in var_re) {
      if(!inherits(data[[var]], "factor")) {
        data[[var]] <- factor(data[[var]])
        warning(paste0("'", var, "' is included as a random effect but is ",
                       "not a factor - changing to factor."))
      }
    }
    
    # Prepare gam() arguments
    gam_args_list <- c(list(formula = update(form, dummy_response ~ .), 
                            data = cbind(dummy_response = 1, data)),
                       gam_args)
    
    # Create matrices based on this formula
    gam_setup <- do.call(what = gam,
                         args = c(gam_args_list, list(fit = FALSE)))
    # Extract column names for design matrices
    term_names <- gam_setup$term.names
    if(is.null(new_data)) {
      Xmat <- gam_setup$X
    } else {
      # Get design matrix for new data set. predict.gam() with
      # type = "lpmatrix" uses none of the quantities that fitting produces,
      # so the unfitted setup can be wrapped in a shell and used instead of
      # fitting a gam to the dummy response, as this used to do.
      Xmat <- predict(gam_shell(gam_setup), newdata = new_data,
                      type = "lpmatrix")
    }
    
    # Fixed effects design matrix
    X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
    subnames_fe <- paste0(names[k], ".", term_names[1:gam_setup$nsdf])
    names_fe <- c(names_fe, subnames_fe)
    
    # Random effects design matrix
    X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
    if(ncol(X_list_re[[k]]) > 0) {
      subnames_re <- paste0(names[k], ".", term_names[-(1:gam_setup$nsdf)])
      names_re <- c(names_re, subnames_re)                    
    }
    
    # Smoothing matrix
    S_list[[k]] <- bdiag_check(gam_setup$S)
    
    # One generalised determinant per penalty matrix, rather than one for the
    # block diagonal of all of this formula's penalties. The likelihood adds
    # -0.5 * log|S_i|+ for each penalty separately, so a linear predictor with
    # several smooths needs them apart; the block diagonal gave their sum to
    # the first smooth and nothing to the rest.
    # vapply, not sapply: a formula with no smooth has an empty penalty list,
    # and sapply() would return a list and coerce log_det_S along with it
    log_det_S <- c(log_det_S, vapply(gam_setup$S, gdeterminant, numeric(1)))
    
    # Number of columns for fixed effects
    ncol_fe <- c(ncol_fe, gam_setup$nsdf)
    
    if(length(gam_setup$smooth) > 0) {
      sub_ncol_re <- matrix(1, nrow = 2, ncol = length(gam_setup$S))
      colnames(sub_ncol_re) <- 1:ncol(sub_ncol_re)
      start_s <- 1
      for (s in 1:length(gam_setup$smooth)) {
        sm <- gam_setup$smooth[[s]]
        # how many penalties for this smooth?
        npen <- length(sm$S)
        # how many parameters for this smooth? 
        npar <- ncol(sm$S[[1]])
        # where does this smooth's parameters start and end?
        sub_ncol_re[, (start_s:(start_s + npen - 1))] <- c(start, start + npar - 1)
        colnames(sub_ncol_re)[start_s:(start_s + npen - 1)] <- rep(sm$label, npen)
        # get names of smooth terms
        # regex from datascience.stackexchange.com/questions/8922
        s_terms <- gsub("(.*)\\..*", "\\1", names_re[sub_ncol_re[1, s]:sub_ncol_re[2, s]])
        s_label <- unique(s_terms)
        names_ncol_re <- c(names_ncol_re, rep(s_label, npen))
        
        # mgcv's L convention: a smooth may combine several penalties through
        # fewer smoothing parameters, with log(lambda) = L * theta. L is the
        # identity for an ordinary smooth, which has one of each.
        L_sm <- if(is.null(sm$L)) diag(npen) else as.matrix(sm$L)
        L_list <- c(L_list, list(L_sm))
        ntheta <- ncol(L_sm)
        sp_names <- c(sp_names, if(is.null(sm$theta.names)) rep(s_label, ntheta)
                      else paste0(s_label, ".", sm$theta.names))
        # A smooth may supply its own starting values, on the log scale.
        # An ordinary smooth supplies none, which means 0, i.e. lambda = 1.
        theta_start <- c(theta_start,
                         if(is.null(sm$theta.start)) rep(0, ntheta)
                         else rep(sm$theta.start, length = ntheta))
        
        start <- start + npar
        start_s <- start_s + npen
      }
      ncol_re <- cbind(ncol_re, sub_ncol_re)
    }    
  }
  colnames(ncol_re) <- names_ncol_re
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  colnames(X_fe) <- names_fe
  X_re <- bdiag_check(X_list_re)
  colnames(X_re) <- names_re
  S <- bdiag_check(S_list)
  L <- bdiag_check(L_list)
  if(!is.null(L)) L <- as.matrix(L)
  
  return(list(X_fe = X_fe, 
              X_re = X_re, 
              S = S,
              log_det_S = log_det_S,
              X_list_fe = X_list_fe, 
              X_list_re = X_list_re, 
              S_list = S_list, 
              ncol_fe = ncol_fe, 
              ncol_re = ncol_re,
              L = L,
              # sp_names stays NULL when there are no smooths, as ncol_re is,
              # so that update_lambda() leaves an empty lambda matrix unnamed.
              # theta_start is numeric(0) instead, so that exp() accepts it.
              sp_names = sp_names,
              theta_start = as.numeric(theta_start)))
}

#' Shell gam object for building prediction matrices
#' 
#' \code{mgcv::predict.gam()} with \code{type = "lpmatrix"} uses only the model
#' frame, the terms objects, the factor levels and contrasts, and the smooth
#' objects -- all of which \code{mgcv::gam(fit = FALSE)} already produces. This
#' wraps that unfitted setup in something \code{predict.gam()} accepts, so that
#' a prediction matrix can be built without fitting anything.
#' 
#' The coefficients are set to zero. \code{predict.gam()} requires the slot to
#' exist and to have the right length, but for \code{type = "lpmatrix"} it
#' returns the design matrix itself and never multiplies by them.
#' 
#' @param G Output of \code{mgcv::gam()} called with \code{fit = FALSE}
#' 
#' @return An object of class "gam", only usable for
#' \code{predict(type = "lpmatrix")}
gam_shell <- function(G) {
  shell <- G[c("pterms", "terms", "smooth", "nsdf", "assign", "xlevels", 
               "contrasts", "pred.formula")]
  shell$model <- G$mf
  shell$na.action <- attr(G$mf, "na.action")
  shell$coefficients <- stats::setNames(rep(0, ncol(G$X)), G$term.names)
  class(shell) <- c("gam", "glm", "lm")
  return(shell)
}
