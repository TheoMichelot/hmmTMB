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
#'   \item obsnames: columns names of data streams that are observations
#'   \item dists: distributions for each data stream
#' }
#'
#' Methods include:
#' \itemize{
#'  \item  NA
#' }
#'

Observation <- R6Class("Observation",

  public = list(
    initialize = function(data, obsnames, dists, par) {
       private$data_ <- data
       private$obsnames_ <- obsnames
       private$dists_ <- dists
       private$par_ <- par
       private$tpar_ <- log(par)
    },

    # Accessors
    data = function() {return(private$data_)},
    obsnames = function() {return(private$obsnames_)},
    dists = function() {return(private$dists_)},
    par = function() {return(private$par_)},
    tpar = function() {return(private$tpar_)}

  ),

  private = list(
    data_ = NULL,
    obsnames_ = NULL,
    dists_ = NULL,
    par_ = NULL,
    tpar_ = NULL

  )
)
