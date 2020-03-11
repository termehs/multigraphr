#' @title Ordered n-tuples of non-negative integers summing to k
#' @description Finds ordered n-tuples of non-integers summing to k.
#' Only practical for n < 15
#' @param n A positive integer
#' @param k A positive integer
#' @return  A matrix X with nchoosek(k+n-1,n-1) rows and n columns. Each row comprises non-negative integers summing to k.
#' @details  To be completed
#' @author Termeh Shafie
#' @examples
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
