
#' R6 class for HMM data
#'
#' Contains a data frame containing observations for the response variables,
#' as well as the covariates.
HMMData <- R6Class(
  classname = "HMMData",
  
  public = list(
    #################
    ## Constructor ##
    #################
    #' @description Create new HMMData object
    #' 
    #' @param data A data frame with at least one column of observations. A column
    #' called "time" can be given for the time of each observation, otherwise
    #' it is assumed that observations are regular in time. Other columns should
    #' include "ID" (time series ID) if necessary, and any covariates needed in
    #' the model.
    #' @param interval Time interval of observation. Optional.
    #' 
    #' @return A new HMMData object
    initialize = function(data, interval = NA) {
      # Check arguments
      private$check_args(data = data, interval = interval)
      
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
    #' @description Data frame
    data = function() {return(private$data_)},
    
    ###################
    ## Other methods ##
    ###################
    #' @description Get IDs
    #' 
    #' @return Vector of time series IDs. Defaults to a vector
    #' of 1s if there is no ID columns in the input data set.
    ID = function() {
      if(is.null(self$data()$ID)) {
        return(factor(rep(1, nrow(self$data()))))
      } else {
        return(self$data()$ID)
      }
    }
  ),
  
  private = list(
    ################
    ## Attributes ##
    ################
    data_ = NULL,
    
    #################################
    ## Check constructor arguments ##
    #################################
    # (For argument description, see constructor)
    check_args = function(data, interval) {
      if(!is.data.frame(data)) {
        stop("'data' must be a data frame")
      }
    }
  )
)
