
#' logLik function for SDE objects
#' 
#' This function makes it possible to call generic R methods such
#' as AIC and BIC on HMM objects. It is based on the number of
#' degrees of freedom of the *conditional* AIC (rather than
#' marginal AIC), i.e., including degrees of freedom from the
#' smooth/random effect components of the model.
#' 
#' @param object HMM model object
#' @param ... For compatibility with S3 method
#' 
#' @return Maximum log-likelihood value for the model, with attributes
#' \code{df} (degrees of freedom) and \code{nobs} (number of observations)
#' 
#' @export
logLik.HMM <- function(object, ...) {
  if(nrow(object$obs()$coeff_re()) + nrow(object$hid()$coeff_re()) > 0) {
    warning("AIC/BIC functions are experimental for models with random effects",
            " or splines. Use at your own risk.")
  }
  
  par_all <- c(object$tmb_rep()$par.fixed, 
               object$tmb_rep()$par.random)
  val <- - object$tmb_obj_joint()$fn(par_all)
  attributes(val)$df <- object$edf()
  attributes(val)$nobs <- nrow(object$obs()$data())
  class(val) <- "logLik"
  return(val)
}
