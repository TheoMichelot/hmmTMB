% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLik.R
\name{logLik.HMM}
\alias{logLik.HMM}
\title{logLik function for SDE objects}
\usage{
\method{logLik}{HMM}(object, ...)
}
\arguments{
\item{object}{HMM model object}

\item{...}{For compatibility with S3 method}
}
\value{
Maximum log-likelihood value for the model, with attributes
\code{df} (degrees of freedom) and \code{nobs} (number of observations)
}
\description{
This function makes it possible to call generic R methods such
as AIC and BIC on HMM objects. It is based on the number of
degrees of freedom of the *conditional* AIC (rather than
marginal AIC), i.e., including degrees of freedom from the
smooth/random effect components of the model.
}
