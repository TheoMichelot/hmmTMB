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
# file description: data class
#
################################################################################

#' Hidden Markov model data class
#'
#' @description Encapsulates data observing one or more time series.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item data: a data frame with at least one column of observations. A column called "time" can be
#'   given for the time of each observation, otherwise assumed observations are regular in time.
#' }
#'
#' Methods include:
#' \itemize{
#'  \item
#' }
#'

HmmData <- R6Class("HmmData",

  public = list(
    initialize = function(data, interval = NA) {
      
      # If time column and interval provided, insert NAs to obtain regular time grid
      if(!is.na(interval) & !is.null(data$time)) {
        # Initialise regularised data set
        data_reg <- NULL
        
        # Loop over IDs
        for(id in unique(data$ID)) {
          sub_data <- subset(data, ID == id)
          sub_n <- nrow(sub_data)
          
          # Regular time grid
          reg_times <- data.frame(time = seq(from = sub_data$time[1], 
                                             to = sub_data$time[sub_n], 
                                             by = interval))
          
          # Insert NAs
          sub_data_reg <- merge(sub_data, reg_times, by = "time", all = TRUE)
          
          # Remove NAs in ID column
          sub_data_reg$ID <- id
          
          # Append to regularised data set
          data_reg <- rbind(data_reg, sub_data_reg)
        }
        
        private$data_ <- data_reg        
      } else {
        private$data_ <- data
      }
    },

    data = function() {return(private$data_)}
  ),

  private = list(
    data_ = NULL
  )
)
