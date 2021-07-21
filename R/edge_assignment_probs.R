#' @title Edge assignment probabilities under the IEAS and ISA model
#' @description Calculates the edge assignment probabilities given a degree sequence under the two ways in which the RSM
#' model can be approximated by the IEA model. This is done by either the IEAS (independent edge assignment of stubs) or
#' the ISA (independent stub assignment) model.
#' @param m integer giving number of edges in multigraph
#' @param deg.seq vector of integers with sum 2m representing the degree sequence of the multigraph
#' @param model character string, either IEAS or ISA
#' @return  A numeric vector Q.seq representing the edge assignment probabilities
#' to all possible vertex pair sites
#' @details The IEAS and ISA edge assignment probabilities to
#' all possible vertex pairs are calculated given a fixed degree sequence \emph{deq.seq}
#' under the IEAS model, and \emph{deg.seq/2m} under the ISA model.
#' Number of possible vertex pair sites to which edges are assigned is given by
#' \eqn{(n+1)n/2} where \emph{n} is number of vertices.
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' ## Under the IEAS model with 10 possible vertex pair sites (4 vertices)
#' edge_assignment_probs(m = 8, deg.seq = c(4,4,4,4), model = 'IEAS')
#'
#' #' ## Under the ISA model with 21 possible vertex pair sites (6 vertices)
#' edge_assignment_probs(m = 10, deg.seq = c(8,4,2,2,2,2), model = 'ISA')
#' @export
#'
edge_assignment_probs <- function(m, deg.seq, model) {
  n <- length(deg.seq)
  if (model == 'IEAS') {
    Q.mat <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i == j) {
          Q.mat[i, j] <- deg.seq[i] * (deg.seq[i] - 1) / (2 * m * (2 * m - 1))
        } else {
          Q.mat[i, j] <- 2 * deg.seq[i] * (deg.seq[j]) / (2 * m * (2 * m - 1))
        }
      }
    }
    Q.seq <-  t(Q.mat)[lower.tri(Q.mat, diag = TRUE)]
    return(Q.seq)
  }
  else if (model == 'ISA') {
    if (sum(deg.seq) == 2 * m) {
      p.seq <- deg.seq / (2 * m)
    } else {
      p.seq <- deg.seq
    }
    Q.mat <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i == j) {
          Q.mat[i, j] <- p.seq[i] ^ 2
        } else {
          Q.mat[i, j] <- 2 * p.seq[i] * p.seq[j]
        }
      }
    }
    Q.seq <-  t(Q.mat)[lower.tri(Q.mat, diag = TRUE)]
    return(Q.seq)
  }
}
