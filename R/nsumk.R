#' @title Ordered \emph{n}-tuples of non-negative integers summing to \emph{k}
#' @description Finds ordered \emph{n}-tuples of non-integers summing to \emph{k}.
#' Only practical for \emph{n} < 15.
#' @param n a positive integer.
#' @param k a positive integer.
#' @return  a matrix with \code{choose(k+n-1,n-1)} rows and \code{n} columns.
#' Each row comprises non-negative integers summing to \code{k}.
#' @details  Useful for finding all possible degree sequences for a network with \code{n} nodes
#' and \code{k}/2 number of edges, or for finding all possible edge multiplicity sequences \code{n} that
#' sum up to \emph{k/2} number of edges.
#' @author Termeh Shafie
#' @seealso \code{\link{gof_sim}}
#' @examples
#' # All possible degree sequences for
#' # a network with 4 nodes and 5 edges
#' D <- nsumk(4, 10)
#'
#' # Remove isolated nodes
#' D <- D[-which(rowSums(D == 0) > 0), ]
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
    warning('n and k are not scalars')
  }
  return(X)
}
