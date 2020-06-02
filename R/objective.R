#' \code{cos.sim.un} calculates the similarity penalty
#'
#' @param Mat is the \eqn{(p \times M)} model matrix
#' @export
cos.sim.un <- function(Mat)
{
  nMat <- t(abs(Mat))
  sim <- nMat# / sqrt(rowSums(nMat * nMat))
  sim <- sim %*% t(sim)
  return(sum(sim[upper.tri(sim)]))
}

#' \code{cos.sim.all} calculates the model similarity matrix
#'
#' @inheritParams cos.sim.un
#' @export
cos.sim.all <- function(Mat)
{
  nMat <- t(abs(Mat))
  sim <- nMat / sqrt(rowSums(nMat * nMat))
  sim <- sim %*% t(sim)
  return(sim)
}

#' \code{cos.sim.max} calculates the maximum model similarity
#'
#' @inheritParams cos.sim.un
#' @export
cos.sim.max <- function(Mat){
  nMat <- t(abs(Mat))
  sim <- nMat / sqrt(rowSums(nMat * nMat))
  sim <- sim %*% t(sim)
  return(max(sim[upper.tri(sim)],na.rm=T))
}
