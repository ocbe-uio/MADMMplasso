#' @title Fit the hierarchical tree structure
#' @description Fit the hierarchical tree structure
#' @param y  N by D matrix of response variables
#' @param h is the tree cut off
#' @return  A trained  tree with the following components:
#' Tree: the tree matrix stored in 1s and 0s
#'  Tw: tree weights associated with the tree matrix. Each weight corresponds to a row in the tree matrix.
#'  h_clust: Summary of the hclust call
#'  y.colnames: names of the response

#' @export
tree_parms <- function(y = y, h = 0.7) {
  m <- ncol(y)
  myDist0 <- 1 - abs(fast_corr(y))
  myDist <- myDist0[lower.tri(myDist0)]
  a0 <- dist(t(y))
  a0[seq_along(a0)] <- myDist
  # hierarchical clustering for multivariate responses
  myCluster_0 <- hclust(a0, method = "complete")
  myCluster <- cbind(ifelse(myCluster_0$merge < 0, -myCluster_0$merge, myCluster_0$merge + m), myCluster_0$height)

  conv0 <- conv_H2T(myCluster, h)
  Tree <- conv0$Tree
  if (is.null(dim(Tree))) {
    Tree <- matrix(Tree, nrow = 1)
  }
  Tw <- conv0$Tw
  idx <- c(apply(Tree, 1, sum) == 1)
  Tree <- Tree[!idx, ]
  if (is.null(dim(Tree))) {
    Tree <- matrix(Tree, nrow = 1)
  }
  Tw <- Tw[!idx]

  out <- list(Tree = Tree, Tw = Tw, h_clust = myCluster_0, y.colnames = colnames(y))

  return(out)
}
