
##' Create block diagonal matrix (safe version)
##' 
##' This version of bdiag checks whether the matrices passed as
##' arguments are NULL. This avoids errors that would arise if
##' using bdiag directly.
##' 
##' @param ... Matrix or list of matrices (only the first argument is used)
##' 
##' @return Block diagonal matrix
##' 
##' @importFrom Matrix bdiag
bdiag_check <- function(...) {
  # Only use first argument (matrix or list of matrices)
  args <- list(...)[[1]]
  
  # Check which matrices are non-NULL
  # (Needs two conditions to keep non-empty vectors which don't have
  # dimensions, and also keep empty matrices with set dimensions)
  check <- sapply(args, function(arg) {
    !is.null(dim(arg)) | length(arg) > 0
  })
  
  if(length(check) == 0)
    return(NULL)
  else
    return(bdiag(args[check]))
}
