#' @title Compute predicted values from a fitted MADMMplasso  object.
#'  Make predictions from a fitted MADMMplasso model
#' @description Compute predicted values from a MADMMplasso  object.
#'  Make predictions from a fitted MADMMplasso model
#' @param object object returned from a call to MADMMplasso
#' @param X  N by p matrix of predictors
#' @param Z  N by nz matrix of modifying variables. These may be observed or
#' the predictions from a supervised learning algorithm that predicts z from
#' test features x  and possibly other features.
#' @param y N by D matrix  of responses.
#' @param lambda  values of lambda at which predictions are desired. If NULL (default), the path of lambda values from the fitted model. are used. If lambda is not NULL, the predictions are made at the closest values to lambda in the lambda path from the fitted model
#' @param ... additional arguments to the generic \code{predict()} method
#' @return  predicted values
#' @export
predict.MADMMplasso <- function(object, X, Z, y, lambda = NULL, ...) {
  lambda.arg <- lambda
  if (is.null(lambda.arg)) {
    lambda <- object$Lambdas[, 1]
    isel <- seq_along(lambda)
  }

  if (!is.null(lambda.arg)) {
    isel <- as.numeric(knn1(matrix(object$Lambdas[, 1], ncol = 1), matrix(lambda.arg, ncol = 1), seq_along(object$Lambdas[, 1])))
  }

  N <- nrow(X)

  p <- ncol(X)
  K <- ncol(as.matrix(Z))
  D <- ncol(y)
  my_W_hat <- generate_my_w(X = X, Z = Z)

  yh <- array(0, c(N, D, length(isel)))
  DEV <- matrix(NA, length(isel))
  pBETA0 <- lapply(
    seq_len(length(isel)),
    function(j) (matrix(0, nrow = (D)))
  )
  pBETA <- lapply(
    seq_len(length(isel)),
    function(j) (matrix(0, nrow = p, ncol = (D)))
  )
  pBETA_hat <- lapply(
    seq_len(length(isel)),
    function(j) (matrix(0, nrow = p + p * K, ncol = (D)))
  )

  pTHETA0 <- lapply(
    seq_len(length(isel)),
    function(j) (matrix(0, nrow = K, ncol = (D)))
  )

  pTHETA <- lapply(
    seq_len(length(isel)),
    function(j) (array(0, c(p, K, (D))))
  )

  ii <- 0
  for (m in isel) {
    ii <- ii + 1
    z <- m
    n_i <- lapply(
      seq_len(max(1)),
      function(j) (matrix(0, nrow = N))
    )
    beta0 <- object$beta0[[z]]
    beta <- as.matrix(object$beta[[z]])
    beta_hat <- as.matrix(object$BETA_hat[[z]])
    theta <- as.array(object$theta[[z]])

    theta0 <- object$theta0[[z]]

    pBETA0[[ii]] <- beta0
    pBETA[[ii]] <- beta
    pBETA_hat[[ii]] <- beta_hat
    pTHETA[[ii]] <- theta
    pTHETA0[[ii]] <- theta0


    n_i <- (model_p(beta0, theta0, beta = beta_hat, X = my_W_hat, Z))
    Dev <- (sum(as.vector((y - (n_i))^2))) / (D * N)
    DEV[ii] <- (Dev)
    yh[, , ii] <- as.matrix(n_i)
  }
  out <- list(y_hat = yh, beta0 = pBETA0, beta = pBETA, beta_hat = pBETA_hat, theta0 = pTHETA0, theta = pTHETA, deviance = DEV)
  return(out)
}
