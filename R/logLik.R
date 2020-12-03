#' @description logLik function for HMM TMB 
#' @param object the HMM model object 
#' @export logLik.HMM
logLik.HMM <- function(object, ...) {
  val <- -object$out()$value 
  attributes(val)$df <- object$edf()
  attributes(val)$nobs <- nrow(object$obs()$data())
  class(val) <- "logLik"
  return(val)
}
