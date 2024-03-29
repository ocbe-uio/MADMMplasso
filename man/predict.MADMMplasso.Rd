% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.MADMMplasso.R
\name{predict.MADMMplasso}
\alias{predict.MADMMplasso}
\title{Compute predicted values from a fitted MADMMplasso  object.
Make predictions from a fitted MADMMplasso model}
\usage{
\method{predict}{MADMMplasso}(object, X, Z, y, lambda = NULL, ...)
}
\arguments{
\item{object}{object returned from a call to MADMMplasso}

\item{X}{N by p matrix of predictors}

\item{Z}{N by nz matrix of modifying variables. These may be observed or
the predictions from a supervised learning algorithm that predicts z from
test features x  and possibly other features.}

\item{y}{N by D matrix  of responses.}

\item{lambda}{values of lambda at which predictions are desired. If NULL (default), the path of lambda values from the fitted model. are used. If lambda is not NULL, the predictions are made at the closest values to lambda in the lambda path from the fitted model}

\item{...}{additional arguments to the generic \code{predict()} method
@return  predicted values}
}
\description{
Compute predicted values from a MADMMplasso  object.
Make predictions from a fitted MADMMplasso model
}
