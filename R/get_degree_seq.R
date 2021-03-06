#' @title Degree sequence of a multigraph
#' @description Gives the degree sequence of the adjacency matrix of an observed graph or multigraph
#' @param adj Matrix of integers
#' @param type Equals 'graph' if adjacency matrix is for graphs (default),
#' equals 'multigraph' if it is the equivalence of the adjacency matrix for multigraphs
#' (with the matrix diagonal representing loops being double counted).
#' @return Vector of integers representing the degree sequence
#' @details Gives the degree sequence of the adjacency matrix of an observed graph or multigraph
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#'  ## Adjacency matrix for undirected network with 3 nodes
#'  A <-  matrix(c(0, 1, 2, 1, 2, 1, 2, 1, 2), nrow=3, ncol=3)
#'  # If A represents a graph
#'  get_degree_seq(adj = A, type = 'graph')
#'  # If A represents a multigraph
#'  get_degree_seq(adj = A, type = 'multigraph')
#' @export
#'
get_degree_seq <- function(adj, type = 'multigraph') {
  if (type == 'multigraph') {
    if (sum(diag(adj)) %% 2 == 1)
      stop("not an adjacency matrix for multigraphs
           with diagonal elements double counted,
           consider type 'graph' instead.")
    if (sum(adj) %% 2 == 1)
      stop("sum of adjacency matrix must be even")
    deg.seq <- (rowSums(adj))
  } else if (type == 'graph') {
    deg.seq <- (rowSums(adj) + diag(adj))
  }
  else {
    stop("type must be defined as either 'graph' or 'multigraph'")
  }
  return(deg.seq)
}
