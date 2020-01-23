
# Copyright (c) 2019 Richard Glennie, Théo Michelot
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################
# hmmtmb project: hidden Markov models using template model builder
#
# file description: Observation distribution class
#
################################################################################

#' Distribution class
#'
#' @description Distribution for observed variables.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item name Name of distribution
#'   \item pdf Probability density/mass function of the distribution
#'   \item link Named list of link functions for distribution parameters
#'   \item invlink Named list of inverse link functions for distribution
#'   parameters
#' }
#'
#' Methods include:
#' \itemize{
#'  \item n2w Transforms parameters from natural to working scale
#'  \item w2n Transforms parameters from working to natural scale
#' }
#'

Dist <- R6Class(
  classname = "Dist",
  
  public = list(
    initialize = function(name, pdf, link, invlink) {
      private$name_ <- name
      private$pdf_ <- pdf
      private$link_ <- link
      private$invlink_ <- invlink
    },
    
    # Accessors
    name = function() {return(private$name_)},
    pdf = function() {return(private$pdf_)},
    link = function() {return(private$link_)},
    invlink = function() {return(private$invlink_)},
    
    # Transform parameters from natural to working scale
    n2w = function(par) {
      # Apply link functions to natural parameters
      wpar_list <- Map(function(fn, arg) {fn(arg)}, private$link_, par)
      wpar <- unlist(wpar_list)
      return(wpar)
    },
    
    # Transform parameters from working to natural scale
    w2n = function(wpar) {
      par <- list()
      
      # Number of parameters for this distribution
      n_par <- length(private$link_)
      # Number of states
      n_state <- length(wpar)/n_par
      
      # Loop over parameters and apply inverse link functions
      for(i in 1:n_par) {
        i0 <- (i-1)*n_state + 1
        i1 <- i*n_state
        sub_wpar <- wpar[i0:i1]
        par[[i]] <- private$invlink_[[i]](sub_wpar)
      }
      
      names(par) <- names(private$link_)
      return(par)
    }
  ),
  
  private = list(
    name_ = NULL,
    pdf_ = NULL,
    link_ = NULL,
    invlink_ = NULL
  )
)
