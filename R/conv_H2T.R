conv_H2T <- function(H, w_max) {
  K <- nrow(H) + 1
  Nd <- cbind(rep((K + 1):(2 * K - 1), each = 2), as.vector(t(H[, 1:2])))
  W_norm <- H[, 3] / max(H[, 3])
  conv_Nd2T(Nd, W_norm, w_max)
}
