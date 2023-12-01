#' @export
plot.cv_MADMMplasso <- function(x, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  plot.args <- list(
    x = log(as.matrix(cvobj$lambda[, 1])), y = as.matrix(cvobj$cvm),
    ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = "Error",
    type = "n"
  )

  new.args <- list(...)
  if (length(new.args)) {
    plot.args[names(new.args)] <- new.args
  }
  do.call("plot", plot.args)
  error_bars(log(cvobj$lambda[, 1]), cvobj$cvup, cvobj$cvlo,
    width = 0.01, col = "darkgrey"
  )
  points(log(as.matrix(cvobj$lambda[, 1])), as.matrix(cvobj$cvm),
    pch = 20,
    col = "red"
  )
  axis(
    side = 3, at = log(as.matrix(cvobj$lambda[, 1])), labels = paste(as.matrix(cvobj$nz)),
    tick = FALSE, line = 0
  )

  abline(v = log(cvobj$lambda.min), lty = 3)
  abline(v = log(cvobj$lambda.1se), lty = 3)
  invisible()
}
