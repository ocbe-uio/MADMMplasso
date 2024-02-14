#' @name MADMMplasso
#' @importFrom graphics abline axis lines matplot points segments text
#' @importFrom methods as
#' @importFrom stats dist hclust lm rbinom var
#' @importFrom Matrix sparseMatrix Diagonal bdiag
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParRegistered %dopar% foreach
#' @importFrom spatstat.sparse as.sparse3Darray
#' @importFrom class knn1
#' @importFrom MASS mvrnorm
#' @importFrom Rcpp sourceCpp
#' @useDynLib MADMMplasso, .registration = TRUE
"_PACKAGE"
