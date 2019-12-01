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
# file description: hidden Markov model class
#
################################################################################

#' Hidden Markov model class
#'
#' @description Encapsulates the hidden Markov model.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item obs: a Observations object
#'   \item hidden: a MarkovChain object
#' }
#'
#' Methods include:
#' \itemize{
#'  \item  NA
#' }
#'

Hmm <- R6Class("Hmm",

  public = list(
    initialize = function(obs, hidden) {
       private$obs_ <- obs
       private$hidden_ <- hidden
    },

    # Accessors
    obs = function() {return(private$obs_)},
    hidden = function() {return(private$hidden_)},
    res = function() {
      if (is.null(private$fit_)) stop("Fit model first")
      return(private$fit_)
    },

    # Fitting
    fit = function() {

      tmb_dat <- list(data = self$obs()$data()$data()[,1],
                            n_states = self$hidden()$nstates())

      tmb_par <- list(ltpm = self$hidden()$par(),
                      distpar = self$obs()$tpar())

      obj <- MakeADFun(tmb_dat, tmb_par, dll = "HmmTmb")

      private$fit_ <- do.call(optim, obj)

    }

  ),

  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    fit_ = NULL

  )
)





