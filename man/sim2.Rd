% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MADAMM11.R
\name{sim2}
\alias{sim2}
\title{Simulate data for the model}
\usage{
sim2(p = 500, n = 100, m = 24, nz = 4, rho = 0.4, B.elem = 0.5)
}
\arguments{
\item{p:}{column for X which is the main effect}

\item{n:}{number of observations
#' @param m: number of responses}
}
\value{
The simulated data with the following components:
Beta: matrix of actual beta coefficients  p by D
Theta: a D by p by K array of actual theta coefficients
Y: a N by D matrix of response variables
X: a N by p matrix of covariates
Z: a N by K matrix of modifiers
}
\description{
Simulate data for the model
}
