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
# file description: hidden Markov observation process class
#
################################################################################

#' Hidden Markov observation class
#'
#' @description Encapsulates the observation model from a hidden Markov model.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item data: a HmmData object
#'   \item dists: named list of distributions for each data stream
#'   \item par: list of observation parameters
#' }
#'
#' Methods include:
#' \itemize{
#'  \item plot_dist Plot histogram of observations, overlaid with the pdf
#'  of the specified distribution for that data stream. Helpful to select
#'  initial parameter values for the model fitting.
#' }
#'

Observation <- R6Class("Observation",

  public = list(
    initialize = function(data, dists, par) {
       private$data_ <- data
       private$dists_ <- dists
       private$par_ <- par
       # private$tpar_ <- log(par)
    },

    # Accessors
    data = function() {return(private$data_)},
    dists = function() {return(private$dists_)},
    par = function() {return(private$par_)},
    tpar = function() {return(private$tpar_)},

    # Histogram of observations with overlaid pdf
    plot_dist = function(name, par = NULL) {
      # Extract observed values for relevant variable
      obs <- private$data_$data()[[name]]
      
      # Histogram of observations
      hist(obs, col = "lightgrey", border = "white", prob = TRUE, 
           main = "", xlab = name)
      
      # Create list of arguments for pdf
      grid <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = 1e3)
      args <- list(grid)
      if(!is.null(par)) {
        for(i in 1:length(par))
          args[[i+1]] <- par[i]        
      } else {
        for(i in 1:length(private$par_))
          args[[i+1]] <- private$par_[[i]]
      }
      
      # Add pdf to histogram plot
      points(grid, do.call(private$dists_[[name]]$pdf(), args), type = "l")
    }
  ),

  private = list(
    data_ = NULL,
    obsnames_ = NULL,
    dists_ = NULL,
    par_ = NULL,
    tpar_ = NULL

  )
)

