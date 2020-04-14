#' @title Edge Assignment Probabilities
#' @description Calculates the edge assignment probabilities given a degree sequence under the two ways in which the RSM
#' model can be approxiated by the IEA model. This is done by either the IEAS (independent edge assignment of stubs) or
#' the ISA (independent stub assignment) model.
#' @param deg.seq
#' @return  A vector Q.seq giving the edge assignment probabilities to all possible vertex pair sites
#' @details  To be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#' @examples
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
