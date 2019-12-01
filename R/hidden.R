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
# file description: hidden Markov chain class
#
################################################################################

#' Hidden Markov chain class
#'
#' @description Encapsulates the Markov chain for the hidden component of the HMM.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item structure: does nothing yet. Will be a matrix with an entry of "." on diagonal, a "0" for
#'   transitions that are not allowed, and a formula "~1" for covarites affecting transitions that are to
#'   be estimated.
#'   \item tpm: an initial transition probability matrix.
#' }
#'
#' Methods include:
#' \itemize{
#'  \item  structure: returns the specified structure of the Markov chain
#'  \item tpm: return current transition probability matrix
#'  \item par: return current parameter estimates for transitions
#'  \item nstates: return number of states in Markov chain
#'  \item update_par (newpar): set parameters to newpar
#'  \item update_tpm (newtpm): set transition probability matrix to newtpm
#' }
#'

MarkovChain <- R6Class("MarkovChain",

  public = list(
    initialize = function(structure, tpm) {
      private$structure_ <- structure
      private$nstates_ <- nrow(structure)
      private$tpm_ <- tpm
      private$par_ <- private$tpm2par(tpm)
    },

    # Accessors
    structure = function() {return(private$structure_)},
    tpm = function() {return(private$tpm_)},
    par = function() {return(private$par_)},
    nstates = function() {return(private$nstates_)},

    # Mutators
    update_par = function(newpar) {
      private$par_ <- newpar
      private$tpm_ <- private$par2tpm(newpar)
    },
    update_tpm = function(newtpm) {
      private$tpm_ <- newtpm
      private$par_ <- private$tpm2par(newtpm)
    }

  ),

  private = list(
    structure_ = NULL,
    par_ = NULL,
    tpm_ = NULL,
    nstates_ = NULL,

    check_structure = function() {
      if (!all(diag(private$structure_) == ".")) stop("Diagonal of structure should be '.'")
      return(TRUE)
    },

    tpm2par = function(tpm) {
      ltpm <- log(tpm / diag(tpm))
      par <- ltpm[!diag(self$nstates())]
      return(par)
    },

    par2tpm = function(par) {
      tpm <- diag(self$nstates())
      tpm[!diag(self$nstates())] <- exp(par)
      tpm <- tpm / rowSums(tpm)
      return(tpm)
    }

  )
)





