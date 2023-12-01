#' Compute predicted values from a fitted pliable  object
#'  Make predictions from a fitted pliable lasso model
#' @param object object returned from a call to pliable
#' @param X  N by p matrix of predictors
#' @param Z  n by nz matrix of modifying variables. These may be observed or
#' the predictions from a supervised learning algorithm that predicts z from
#' test features x  and possibly other features.
#' @param y N by D matrix  of responses.
#' @param lambda  TODO: fill in description
#' @param ... additional arguments to the generic \code{predict()} method
#'  @return  predicted values
#' @export
predict.MADMMplasso <- function(object, X, Z, y, lambda = NULL, ...) {
  lambda.arg <- lambda
  if (is.null(lambda.arg)) {
    lambda <- object$Lambdas[, 1]
    isel <- 1:length(lambda)
  }

  if (!is.null(lambda.arg)) {
    isel <- as.numeric(knn1(matrix(object$Lambdas[, 1], ncol = 1), matrix(lambda.arg, ncol = 1), 1:length(object$Lambdas[, 1])))
  }

  N <- nrow(X)

  p <- ncol(X)
  K <- ncol(as.matrix(Z))
  D <- dim(y)[2]
  my_W_hat <- generate_my_w(X = X, Z = Z, quad = TRUE)

  yh <- array(0, c(N, D, length(isel)))
  DEV <- matrix(NA, length(isel))
  my_theta <- array(0, c(ncol(X), ncol(Z), length(isel)))
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

  pY_HAT <- lapply(
    seq_len(length(isel)),
    function(j) (matrix(0, nrow = N, ncol = (D)))
  )

  ii <- 0
  for (m in isel) {
    ii <- ii + 1
    z <- m
    n_i <- lapply(
      seq_len(max(1)),
      function(j) (matrix(0, nrow = N))
    )
    pr <- lapply(
      seq_len(max(1)),
      function(j) (matrix(0, nrow = N))
    )
    beta0 <- object$beta0[[z]]
    beta <- object$beta[[z]]
    beta_hat <- object$BETA_hat[[z]]
    theta <- object$theta[[z]]

    theta0 <- object$theta0[[z]]

    pBETA0[[ii]] <- beta0
    pBETA[[ii]] <- beta
    pBETA_hat[[ii]] <- beta_hat
    pTHETA[[ii]] <- theta
    pTHETA0[[ii]] <- theta0
  

    n_i <- (model_p(beta0, theta0, beta = beta_hat, theta, X = my_W_hat, Z))
    Dev <- (sum(as.vector((y - as.matrix(n_i))^2))) / (D * N)
    DEV[ii] <- (Dev)
    yh[, , ii] <- as.matrix(n_i)
  }
  out <- list(y_hat = yh, beta0 = pBETA0, beta = pBETA, beta_hat = pBETA_hat, theta0 = pTHETA0, theta = pTHETA, deviance = DEV)
  return(out)
}
