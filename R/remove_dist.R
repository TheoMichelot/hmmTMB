
#' Remove distribution from package's stored distributions
#' 
#' @param name name of distribution to remove 
#' @param compile_cpp if FALSE then package is not recompiled after changes 
#' 
#' @export
remove_dist <- function(name, compile_cpp = TRUE) {
  # load list of stored distributions
  load(paste0(find.package("hmmTMB"), "/distnames.RData"))
  if (! (name %in% distnames$name)) stop("No distribution stored by that name.")
  # get package root 
  root <- find.package("hmmTMB")
  # delete distribution definition file 
  file.remove(paste0(root, "/include/dist__", name, ".hpp"))
  # remove from added dists file 
  added_dist_file <- scan(paste0(root, "/include/added_dists.hpp"), character(), sep = "\n")
  find_entry <- stringr::str_detect(added_dist_file, name)
  added_dist_file <- added_dist_file[!find_entry]
  cat(added_dist_file, file = paste0(root, "/include/added_dists.hpp"), sep = "\n")
  # remove from dist file 
  dist_file <- scan(paste0(root, "/include/dist.hpp"), character(), sep = "\n")
  find_entry <- which(stringr::str_detect(dist_file, name))
  dist_file <- dist_file[-c(find_entry - 1, find_entry)]
  cat(dist_file, file = paste0(root, "/include/dist.hpp"), sep = "\n")
  # remove from distnames 
  distnames <- distnames[distnames$name != name,]
  save(distnames, file = paste0(root, "/distnames.RData"))
  # recompile package
  if (compile_cpp) {
    # create a dummy compilation file 
    cat('#include "hmmTMB.hpp"\n', file = paste0(root, "/include/hmmTMB.cpp"))
    # compile using dummy file 
    comp <- tryCatch(TMB::compile(paste0(root, "/include/hmmTMB.cpp")))
    # copy compiled library to lib folder
    file.copy(from = paste0(root, "/include/hmmTMB.so"), to = paste0(root, "/libs/"), overwrite = TRUE)
    # restart R 
    if (exists(".rs.restartR")){
      .rs.restartR() 
    } else {
      warning("You must restart the R session for the new distribution to be available.")
    }
  }
  invisible(c(name, compile_cpp))
}
