#' @title Degree sequence of multigraph
#' @description gives the degree sequence of observed graph or multigraph
#' @param adj integer matrix
#' @param type equals 'graph' if adjacency matrix is for graphs (default),
#' equals 'multigraph' if it is the equivalence of adjacency matrix for multigraphs
#' (with diagonal double counted)
#' @return vector of integers representing the degree sequence
#' @details  to be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#' @examples
#'

get_degree_seq <- function(adj, type = 'multigraph') {
  if(type == 'multigraph'){
    deg.seq <- sort(rowSums(adj))
  } else if(type == 'graph'){
    deg.seq <- sort(rowSums(adj)+diag(adj))
  }
  return(deg.seq)
}
