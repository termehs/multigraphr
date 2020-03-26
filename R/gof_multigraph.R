#' @title Goodness of Fit Tests for Random Multigraph Models
#' @description gives the degree sequence of observed graph or multigraph
#' @param adj integer matrix
#' @param type equals 'graph' if adjacency matrix is for graphs (default),
#' equals 'multigraph' if it is the equivalence of adjacency matrix for multigraphs
#' (with diagonal double counted)
#' model = either RSM, IEAS or ISA with modelled degree sequence deg.mod
#' vechyp= hypothetical degree sequence,
#' hypothetical degree sequence:
#'   if 'IEA': d = deg.hyp fully specified
#'   if 'ISA': p = deg.hyp/2m fully specified
#'   if 'IEA' and deg.hyp = 0, d estimated from edge mutliplicity sequence
#'   if 'ISA' and deg.hyp = 0, p estimated from edge mutliplicity sequence
#' @return vector of integers representing the degree sequence of the multigraph
#' @details  to be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#' @examples
#'
gof_multigraph <- function(m, model, deg.mod, hyp, deg.hyp) {
  n <- length(deg.mod)
  r <- choose(n + 1, 2)

  # model specification: IEAS or ISA
  if (model == 'IEAS') {
    Q.seq <- edge_assignment_probs(m, deg.mod, model = 'IEAS')
    # use IEA model to find probabilities of each multiplicity sequence/multigraph
    m.seq <- nsumk(r, m)
    prob.mg <-
      apply(m.seq, 1, function(x)
        dmultinom(x, prob = Q.seq))
  }
  else if (model == 'ISA') {
    Q.seq <- edge_assignment_probs(m, deg.mod, model = 'ISA')
    # use ISA model to find probabilities of each multiplicity sequence/multigraph
    m.seq <- nsumk(r, m)
    prob.mg <-
      apply(m.seq, 1, function(x)
        dmultinom(x, prob = Q.seq))
  }
  else if (model == 'RSM') {
    # use RSM model to find probabilities of each multiplicity sequence/multigraph
    rsm <- rsm_model(deg.mod)
    m.seq <- rsm$m.seq
    prob.mg <- rsm$prob.dists[, 1]
  }

  # according to hypothesis, get the values of test statistics
  # hypothesis IEAS
  if (hyp == 'IEAS') {
    if (sum(deg.hyp) > 0) {
      dof <- r - 1
      # fully specified with deg.hyp
      Q.seq <- edge_assignment_probs(m, deg.hyp, model = 'IEAS')
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
    else if (sum(deg.hyp) == 0) {
      dof <- r - n
      # estimated dest from data
      deg.est <- matrix(0, nrow(m.seq), n)
      for (i in 1:nrow(m.seq)) {
        M <- matrix(0, n, n)
        M[lower.tri(M, diag = TRUE)] <- 1
        M[M == 1] <- m.seq[i,]
        M <- M + t(M)
        deg.est[i,] <- colSums(M)
      }
      Q.seq <- matrix(0, nrow(m.seq), r)
      for (d in 1:nrow(deg.est)) {
        deg.seq <- deg.est[d,]
        Q.seq[d,] <-
          edge_assignment_probs(m, deg.seq, model = 'IEAS')
      }
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
  }
  # hypothesis ISA
  else if (hyp == 'ISA') {
    if (sum(deg.hyp) > 0) {
      dof <- r - 1
      # fully specified with deg.hyp
      Q.seq <- edge_assignment_probs(m, deg.hyp, model = 'ISA')
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
    else if (sum(deg.hyp) == 0) {
      dof <- r - n
      # estimated dest from data
      deg.est <- matrix(0, nrow(m.seq), n)
      for (i in 1:nrow(m.seq)) {
        M <- matrix(0, n, n)
        M[lower.tri(M, diag = TRUE)] <- 1
        M[M == 1] <- m.seq[i,]
        M <- M + t(M)
        deg.est[i,] <- colSums(M) / (2 * m)
      }
      Q.seq <- matrix(0, nrow(m.seq), r)
      for (d in 1:nrow(deg.est)) {
        deg.seq <- deg.est[d,]
        Q.seq[d,] <-
          edge_assignment_probs(m, deg.seq, model = 'ISA')
      }
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
  }
}
