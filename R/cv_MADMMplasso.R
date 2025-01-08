#' @title Carries out cross-validation for a MADMMplasso model over a path of regularization values
#' @description Carries out cross-validation for a MADMMplasso model over a path of regularization values
#' @inheritParams MADMMplasso
#' @param fit  object returned by the MADMMplasso function
#' @param nfolds  number of cross-validation folds
#' @param foldid  vector with values in 1:K, indicating folds for K-fold CV. Default NULL
#' @param lambda  user specified lambda_3 values.
#' @param rho the Lagrange variable for the ADMM. This value is updated during the ADMM call based on a certain condition.
#' @param TT The results from the hierarchical clustering of the response matrix.
#' This should same as the parameter tree used during the MADMMplasso call.
#' @return  results containing the CV values
#' @example inst/examples/cv_MADMMplasso_example.R
#' @export
cv_MADMMplasso <- function(fit, nfolds, X, Z, y, alpha = 0.5, lambda = fit$Lambdas, max_it = 50000, e.abs = 1E-3, e.rel = 1E-3, nlambda, rho = 5, my_print = FALSE, alph = 1, foldid = NULL, pal = cl == 1L, gg = c(7, 0.5), TT, tol = 1E-4, cl = 1L, legacy = FALSE) {
  BIG <- 10e9
  no <- nrow(X)
  ggg <- vector("list", nfolds)

  yhat <- array(NA, c(no, ncol(y), length(lambda[, 1])))

  if (is.null(foldid)) {
    foldid <- sample(rep(1:nfolds, ceiling(no / nfolds)), no, replace = FALSE)
  }

  nfolds <- length(table(foldid))

  for (ii in 1:nfolds) {
    if (my_print) cat("fold,", ii)
    oo <- foldid == ii

    ggg[[ii]] <- MADMMplasso(X = X[!oo, , drop = FALSE], Z = Z[!oo, , drop = FALSE], y = y[!oo, , drop = FALSE], alpha = alpha, my_lambda = lambda, lambda_min = 0.01, max_it = max_it, e.abs = e.abs, e.rel = e.rel, nlambda = length(lambda[, 1]), rho = rho, tree = TT, my_print = my_print, alph = alph, pal = pal, gg = gg, tol = tol, cl = cl, legacy)

    cv_p <- predict.MADMMplasso(ggg[[ii]], X = X[oo, , drop = FALSE], Z = Z[oo, ], y = y[oo, ])
    ggg[[ii]] <- 0
    yhat[oo, , seq(nlambda)] <- cv_p$y_hat[, , seq(nlambda)]
  }

  ym <- array(y, dim(yhat))
  err <- apply((ym - yhat)^2, c(1, 3), sum)

  mat <- err
  outmat <- matrix(NA, nfolds, ncol(mat))
  mat[is.infinite(mat)] <- NA
  for (i in seq(nfolds)) {
    mati <- mat[foldid == i, , drop = FALSE]
    outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
  }

  non_zero <- c(fit$path$nzero)

  cvm <- (apply(err, 2, mean, na.rm = TRUE)) / ncol(y)
  nn <- apply(!is.na(err), 2, sum, na.rm = TRUE)
  cvsd <- sqrt(apply(err, 2, var, na.rm = TRUE) / (ncol(y) * nn))

  cvm.nz <- cvm
  cvm.nz[non_zero == 0] <- BIG
  imin <- which.min(cvm.nz)
  imin.1se <- which(cvm < cvm[imin] + cvsd[imin])[1]

  out <- list(
    lambda = fit$Lambdas, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
    cvlo = cvm - cvsd, nz = c(fit$path$nzero),
    lambda.min = fit$Lambdas[imin, 1], lambda.1se = fit$Lambdas[imin.1se, 1]
  )
  class(out) <- "cv_MADMMplasso"

  return(out)
}
