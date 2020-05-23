
#' Hidden Markov model data class
#'
#' @description Encapsulates data observing one or more time series.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item data: a data frame with at least one column of observations. 
#'   A column called "time" can be given for the time of each observation, 
#'   otherwise assumed observations are regular in time. Other columns should
#'   include "ID" (time series ID), and any covariates needed in the model.
#' }
#'
#' @section Methods:
#' \itemize{
#'   
#'   \item{\code{data()}}{data frame}
#'   
#'   \item{\code{ID()}}{vector of time series IDs (vector of 1s if not provided
#'  in input data frame)}
#' }

HmmData <- R6Class(
  classname = "HmmData",
  
  public = list(
    #################
    ## Constructor ##
    #################
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
    
    ###############
    ## Accessors ##
    ###############
    data = function() {return(private$data_)},
    
    ###############################
    ## Vector of time series IDs ##
    ###############################
    ID = function() {
      if(is.null(self$data()$ID)) {
        return(factor(rep(1, nrow(self$data()))))
      } else {
        return(self$data()$ID)
      }
    }
  ),
  
  private = list(
    data_ = NULL
  )
)
