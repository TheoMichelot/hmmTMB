
#' R6 class for probability distribution
#' 
#' Contains the probability density/mass function, and the link and inverse link 
#' functions for a probability distribution.
#' 
#' @export
Dist <- R6Class(
  classname = "Dist",
  
  public = list(

    # Constructor -------------------------------------------------------------
    #' @description Create a Dist object
    #' 
    #' @param name Name of distribution
    #' @param pdf Probability density/mass function of the distribution
    #' (e.g. \code{dnorm} for normal distribution).
    #' @param rng Random generator function of the distribution (e.g.
    #' \code{rnorm} for normal distribution).
    #' @param link Named list of link functions for distribution parameters
    #' @param invlink Named list of inverse link functions for distribution
    #' parameters
    #' @param npar Number of parameters of the distribution
    #' @param parnames character vector with name of each parameter 
    #' @param parapprox Function that takes a sample and produces approximate values for the unknown parameters
    #' @param fixed vector with element for each parameter which is TRUE if parameter is fixed
    #' @param cpp location of C++ TMB definition of distribution
    #' @param compile_cpp if FALSE then distribution code in cpp is added but not compiled into package
    #' 
    #' @return A new Dist object
    initialize = function(name, 
                          pdf, 
                          rng, 
                          link, 
                          invlink, 
                          npar, 
                          parnames, 
                          parapprox = NULL, 
                          fixed = NULL, 
                          cpp = NULL, 
                          compile_cpp = TRUE) {
      # Check arguments
      private$check_args(name = name, 
                         pdf = pdf, 
                         rng = rng, 
                         link = link, 
                         invlink = invlink, 
                         npar = npar, 
                         parnames = parnames)
      
      # Define private data members 
      private$name_ <- name
      private$pdf_ <- pdf
      private$rng_ <- rng
      private$link_ <- link
      private$invlink_ <- invlink
      private$npar_ <- npar
      private$parnames_ <- parnames 
      private$parapprox_ <- parapprox 
      private$fixed_ <- fixed
      
      # All parameters are unfixed by default
      if (is.null(fixed)) private$fixed_ <- rep(FALSE, npar)
      
      # Load list of distributions included in the package and check against 
      # specified distribution name 
      load(paste0(find.package("hmmTMB"), "/distnames.RData"))
      if(!name %in% distnames$name) {
        if (is.null(cpp)) {
          stop(paste0("must provide 'cpp' argument for new distributions or use one of the
                      existing distributions:'", 
                      paste(distnames$name, sep = "", collapse = "', '"), 
                      "'"))
        } else { 
          # add distribution to package 
          warning("Adding new distribution to package. If you want to delete this 
                  distribution later, use the function remove_dist().")
          # get package root 
          root <- find.package("hmmTMB")
          # copy C++ TMB definition to package root 
          file.copy(from = cpp, to = paste0(root, "/include/dist__", name, ".hpp"))
          # include in added_dists file 
          added_dist_file <- scan(paste0(root, "/include/added_dists.hpp"), character(), sep = "\n")
          new_dist_file <- c(added_dist_file[-length(added_dist_file)],
                             paste0("#include \"dist__", name, ".hpp\""), 
                             added_dist_file[length(added_dist_file)])
          cat(new_dist_file, file = paste0(root, "/include/added_dists.hpp"), sep = "\n")
          # add to list of available TMB distributions 
          dist_file <- scan(paste0(root, "/include/dist.hpp"), character(), sep = "\n")
          len <- length(dist_file)
          code <- max(distnames$code) + 1 
          new_dist_file <- c(dist_file[-((len - 4):len)], 
                             paste0("  case ", code, ":"), 
                             paste0("    return(std::unique_ptr<Dist<Type>>(new ", name, "<Type>));"), 
                             dist_file[(len - 4):len])
          cat(new_dist_file, file = paste0(root, "/include/dist.hpp"), sep = "\n")
          
          # recompile TMB library 
          if (compile_cpp) {
            # add to distnames data 
            distnames <- rbind(distnames, c(name, code))
            distnames$code <- as.numeric(distnames$code)
            save(distnames, file = paste0(root, "/distnames.RData"))
            # create a dummy compilation file 
            cat('#include "hmmTMB.hpp"\n', file = paste0(root, "/include/hmmTMB.cpp"))
            # compile using dummy file 
            comp <- tryCatch(TMB::compile(paste0(root, "/include/hmmTMB.cpp"), flags = "-Wno-ignored-attributes"))
            if ("try-error" %in% class(comp)) {
              cat(new_dist_file, file = paste0(root, "/include/added_dists.hpp"), sep = "\n")
              cat(new_dist_file, file = paste0(root, "/include/dist.hpp"), sep = "\n")
            }
            # copy compiled library to lib folder
            file.copy(from = paste0(root, "/include/hmmTMB.so"), to = paste0(root, "/libs/"), overwrite = TRUE)
            # restart R 
            if (exists(".rs.restartR")){
              .rs.restartR() 
            } else {
              warning("You must restart the R session for the new distribution to be available.")
            }
          }
        }
        
      } else {
        # load code unique to this distribution 
        private$code_ <- distnames$code[which(distnames$name == name)]       
      }
    },
    
    # Accessors ----------------------------------------------------------------
    #' @description Return name of Dist object
    name = function() {return(private$name_)},
    
    #' @description Return pdf of Dist object
    pdf = function() {return(private$pdf_)},
    
    #' @description Return random generator function of Dist object
    rng = function() {return(private$rng_)},
    
    #' @description Return link function of Dist object
    link = function() {return(private$link_)},
    
    #' @description Return inverse link function of Dist object
    invlink = function() {return(private$invlink_)},
    
    #' @description Return number of parameters of Dist object
    npar = function() {return(private$npar_)},
    
    #' @description Return names of parameters 
    parnames = function() {return(private$parnames_)}, 
    
    #' @description Return function that approximates parameters
    parapprox = function() {return(private$parapprox_)}, 
  
    #' @description Return which parameters are fixed 
    fixed = function() {return(private$fixed_)}, 
    
    #' @description Return code of Dist object
    code = function() {return(private$code_)},
  
    # Mutators ----------------------------------------------------------------

    #' @description Set number of parameters this distribution has
    #' 
    #' @param new_npar Number of parameters
    set_npar = function(new_npar) {private$npar_ <- new_npar}, 
    
    #' @description Set parameter names 
    #' 
    #' @param new_parnames Parameter names
    set_parnames = function(new_parnames) {private$parnames_ <- new_parnames}, 

    # Other methods -----------------------------------------------------------

    #' @description Evaluate probability density/mass function
    #' 
    #' This method is used in Observation$obs_probs. It is a wrapper around self\$pdf,
    #' which prepares the parameters and passes them to the function.
    #' 
    #' @param x Value at which the function should be evaluated
    #' @param par Vector of parameters. The entries should be named if
    #' they are not in the same order as expected by the R function. (E.g.
    #' shape/scale rather than shape/rate for gamma distribution.)
    #' @param log Logical. If TRUE, the log-density is returned. 
    #' Default: FALSE.
    #' 
    #' @return Probability density/mass function evaluated at x for
    #' parameters par
    pdf_apply = function(x, par, log = FALSE) {
      args <- list(x = x)
      args <- c(args, par, log = log)
      do.call(self$pdf(), args)
    },
    
    #' @description Random number generator
    #' 
    #' This method is a wrapper around self\$rng, which prepares the parameters and
    #' passes them to the function.
    #' 
    #' @param n Number of realisations to generate
    #' @param par Vector of parameters. The entries should be named if
    #' they are not in the same order as expected by the R function. (E.g.
    #' shape/scale rather than shape/rate for gamma distribution.)
    #' 
    #' @return Vector of \code{n} realisations of this distribution
    rng_apply = function(n, par) {
      args <- list(n = n)
      args <- c(args, par)
      do.call(self$rng(), args)
    },
    
    #' @description Natural to working parameter transformation
    #' 
    #' This method transforms parameters from the natural scale (i.e., their domain
    #' of definition) to the "working" or "linear predictor" scale (i.e., the real
    #' line). It is a wrapper for self\$link.
    #' 
    #' @param par List of parameters
    #' 
    #' @return Vector of parameters on the working scale
    n2w = function(par) {
      # Apply list of link functions or an all-in-one link? 
      if (!is.list(self$link())) {
        # Number of states
        n_states <- length(par[[1]])
        wpar <- self$link()(par, n_states)
      } else {
        # Apply link functions to natural parameters
        wpar_list <- Map(function(fn, arg) {fn(arg)}, self$link(), par)
        wpar <- unlist(wpar_list)
      }
      return(wpar)
    },
    
    #' @description Working to natural parameter transformation
    #' 
    #' This method transforms parameters from the "working" or "linear predictor" 
    #' scale (i.e., the real line) to the natural scale (i.e., their domain
    #' of definition). It is a wrapper for self\$invlink.
    #' 
    #' @param wpar Vector of working parameters
    #' @param as_matrix Logical. If TRUE, the natural parameters are
    #' returned as a matrix with one row for each state and one column
    #' for each parameter. If FALSE, the natural parameters are returned
    #' as a list (default).
    #' 
    #' @return List or matrix of parameters on natural scale
    w2n = function(wpar, as_matrix = FALSE) {
      invlink <- self$invlink()
      # Number of parameters for this distribution
      n_par <- self$npar()
      # Number of states
      n_states <- length(wpar)/n_par
      # Apply list of invlink functions or an all-in-one invlink? 
      if (!is.list(invlink)) {
        par_list <- invlink(wpar, n_states)
        nms <- names(par_list)
      } else {
        # Apply inverse link functions
        par_list <- lapply(seq_along(invlink), function(i) {
          ind <- ((i-1) * n_states + 1) : (i * n_states)
          invlink[[i]](wpar[ind])
        })
        nms <- names(invlink)
      }
      
      # Format into matrix or list
      if(as_matrix) {
        par <- matrix(unlist(par_list), nrow = n_states)
        colnames(par) <- names(invlink)        
      } else {
        par <- par_list
        names(par) <- nms        
      }
      
      return(par)
    }
    
  ),
  
  private = list(

    # Data members ------------------------------------------------------------
    name_ = NULL, # name of distribution 
    pdf_ = NULL, # pdf function 
    rng_ = NULL, # random number generator 
    link_ = NULL, # list of link functions or all-in-one link function 
    invlink_ = NULL, # list of inverse link functions or all-in-one inverse link 
    npar_ = NULL, # number of parameters for this distribution 
    parnames_ = NULL,  # name of each parameter 
    parapprox_ = NULL,  # function to compute approximate parameters from sample
    fixed_ = NULL,  # TRUE/FALSE for each parameter on whether it is fixed or estimated
    code_ = NULL, # unique distribution code
  
    # (For argument description, see constructor)
    check_args = function(name, pdf, rng, link, invlink, npar, parnames) {
      if(!is.character(name)) {
        stop("'name' must be a character string")
      }
      
      if(!is.function(pdf)) {
        stop(paste("'pdf' must be an R function for the probability density",
                   "function of the distribution (and not the name of the function",
                   "as a string)."))
      }
      
      if(!is.function(rng)) {
        stop(paste("'rng' must be an R function for the random number generator",
                   "of the distribution (and not the name of the function",
                   "as a string)."))
      }
      
      if(length(link) != length(invlink) | !all(names(link) == names(invlink))) {
        stop(paste("'link' and 'invlink' should be lists of the same length, and with", 
                   "the same names (corresponding to the names of the distribution",
                   "parameters"))
      }
      
      if (is.list(link)) {
        if(!all(sapply(link, is.function)) | !all(sapply(invlink, is.function))) {
          stop(paste("Elements of 'link' and 'invlink' should be R functions"))
        }
      } else {
        if (!is.function(link) | !is.function(invlink)) {
          stop(paste("'link' and 'invlink' should either be functions or lists of 
                     functions"))
        }
      }
      
      if(!is.numeric(npar) | npar < 1) {
        stop("'npar' should be a numeric >= 1")
      }
      
      if (length(parnames) != npar) {
        stop("parnames not same length as number of parameters")
      }
      
    }
  )
)
