#' @title Ordered \emph{n}-tuples of non-negative integers summing to \emph{k}
#' @description Finds ordered \emph{n}-tuples of non-integers summing to \emph{k}.
#' Only practical for \emph{n} < 15.
#' @importFrom utils combn
#' @param n a positive integer.
#' @param k a positive integer.
#' @return  A matrix with \code{choose(k+n-1,n-1)} rows and \code{n} columns.
#' Each row comprises non-negative integers summing to \code{k}.
#' @details  Useful for finding all possible degree sequences for a network with \code{n} nodes
#' and \code{k}/2 number of edges, or for finding all possible
#' edge multiplicity sequence \code{n} that sum up to \emph{k} number of edges.
#' The number of vertex pair sites (or length of edge multiplicity sequence) for
#' a multigraph with \emph{n} nodes is given by \eqn{n(n+1)/2}.
#' @author Termeh Shafie
#' @seealso \code{\link{gof_sim}}
#' @examples
#' # All possible degree sequences for
#' # a network with 4 nodes and 5 edges
#' D <- nsumk(4, 10)
#'
#' # Remove degree sequences with isolated nodes
#' D <- D[-which(rowSums(D == 0) > 0), ]
#'
#' # All edge multiplicity sequences/multigraph with 2 nodes and 4 edges
#' r <- (2*3)/2 # vertex pair sites (or length of edge multiplicity sequences)
#' mg <- nsumk(r,4) # number of rows give number of possible multigraphs
#' @export
#'
nsumk <- function(n, k) {
  if (is.atomic(n) &&
      length(n) == 1L && is.atomic(k) && length(k) == 1L) {
    tot = choose(k + n - 1, n - 1)
    divs = cbind(rep(0, tot), t(combn((1:(
      k + n - 1
    )), n - 1)), cbind(rep(1, tot) * (k + n)))
    X = t(apply(divs, 1, diff)) - 1
  } else{
    stop('n and k must be scalars')
  }
  return(X)
}
