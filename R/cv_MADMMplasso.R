#' Carries out cross-validation for  a  pliable lasso model over a path of regularization values
#' @param fit  object returned by the pliable function
#' @param X  N by p matrix of predictors
#' @param Z N by K matrix of modifying variables. The elements of Z  may
#' represent quantitative or categorical variables, or a mixture of the two.
#'  Categorical varables should be coded by 0-1 dummy variables: for a k-level
#' variable, one can use either k or k-1  dummy variables.
#' @param y N by D-matrix of responses. The X and Z variables are centered in
#' the function. We recommmend that x and z also be standardized before the call
#' @param nfolds  number of cross-validation folds
#' @param foldid  vector with values in 1:K, indicating folds for K-fold CV. Default NULL
#' @param alpha TODO: add parameter description
#' @param lambda TODO: add parameter description
#' @param max_it TODO: add parameter description
#' @param e.abs TODO: add parameter description
#' @param e.rel TODO: add parameter description
#' @param nlambda TODO: add parameter description
#' @param rho TODO: add parameter description
#' @param my_print TODO: add parameter description
#' @param alph TODO: add parameter description
#' @param foldid TODO: add parameter description
#' @param parallel TODO: add parameter description
#' @param pal TODO: add parameter description
#' @param gg TODO: add parameter description
#' @param TT TODO: add parameter description
#' @param tol TODO: add parameter description
#' @param cl TODO: add parameter description
#' @return  results containing the CV values
#' @example inst/examples/cv_MADMMplasso_example.R
#' @export
cv_MADMMplasso <- function(fit, nfolds, X, Z, y, alpha = 0.5, lambda = fit$Lambdas, max_it = 50000, e.abs = 1E-3, e.rel = 1E-3, nlambda, rho = 5, my_print = F, alph = 1, foldid = NULL, parallel = T, pal = 0, gg = c(7, 0.5), TT, tol = 1E-4, cl = 2) {
  BIG <- 10e9
  no <- nrow(X)
  ni <- ncol(X)
  nz <- ncol(Z)
  ggg <- vector("list", nfolds)

  yhat <- array(NA, c(no, dim(y)[2], length(lambda[, 1])))

  my_nzero <- matrix(0, nfolds, length(lambda[, 1]))

  if (is.null(foldid)) {
    foldid <- sample(rep(1:nfolds, ceiling(no / nfolds)), no, replace = FALSE)
  }

  nfolds <- length(table(foldid))

  status.in <- NULL

  for (ii in 1:nfolds) {
    print(c("fold,", ii))
    oo <- foldid == ii

    ggg[[ii]] <- MADMMplasso(X = X[!oo, , drop = F], Z = Z[!oo, , drop = F], y = y[!oo, , drop = F], alpha = alpha, my_lambda = lambda, lambda_min = .01, max_it = max_it, e.abs = e.abs, e.rel = e.rel, nlambda = length(lambda[, 1]), rho = rho, tree = TT, my_print = my_print, alph = alph, cv = T, parallel = parallel, pal = pal, gg = gg, tol = tol, cl = cl)

    cv_p <- predict.MADMMplasso(ggg[[ii]], X = X[oo, , drop = F], Z = Z[oo, ], y = y[oo, ])
    ggg[[ii]] <- 0
    yhat[oo, , seq(nlambda)] <- cv_p$y_hat[, , seq(nlambda)]
  }

  ym <- array(y, dim(yhat))
  err <- apply((ym - yhat)^2, c(1, 3), sum)

  mat <- err
  outmat <- matrix(NA, nfolds, ncol(mat))
  good <- matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] <- NA
  for (i in seq(nfolds)) {
    mati <- mat[foldid == i, , drop = FALSE]
    outmat[i, ] <- apply(mati, 2, mean, na.rm = TRUE)
  }

  non_zero <- c(fit$path$nzero)

  cvm <- (apply(err, 2, mean, na.rm = T)) / dim(y)[2]
  nn <- apply(!is.na(err), 2, sum, na.rm = T)
  cvsd <- sqrt(apply(err, 2, var, na.rm = T) / (dim(y)[2] * nn))

  cvm.nz <- cvm
  cvm.nz[non_zero == 0] <- BIG
  imin <- which.min(cvm.nz)
  imin.1se <- which(cvm < cvm[imin] + cvsd[imin])[1]

  out <- list(lambda = fit$Lambdas, cvm = cvm, cvsd = cvsd, cvup = cvm +
    cvsd, cvlo = cvm - cvsd, nz = c(fit$path$nzero), lambda.min = fit$Lambdas[imin, 1], lambda.1se = fit$Lambdas[imin.1se, 1])
  class(out) <- "cv_MADMMplasso"

  return(out)
}
