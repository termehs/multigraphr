#' @title Goodness of Fit Tests for Random Multigraph Models
#' @description  Goodness of fits tests of multigraph models using Pearson (S) and information divergence (A) tests statistics
#' under the random stub matching (RSM) and by independent edge assignments (IEA) model,
#' where the latter is either independent edge assignments of stubs (IEAS) or independent stub assignment (ISA)
#' @param m number of edges
#' @param model = assumed model, either RSM, IEAS or ISA
#' @param deg.mod = modelled degree sequence, vector of integers with even sum
#' @param hyp = hypothesis, either IEAS or ISA
#' @param deg.hyp = hypothetical degree sequence, vector of integers with even sum
#'  hypothetical degree sequence:
#'   if 'IEAS': simple IEAS hypothesis with fully specified degree sequence deg.hyp
#'   if 'ISA': simple ISA hypothesis with with fully specified stub assignment probabilities deg.hyp/2m
#'   if 'IEAS': and deg.hyp = 0: composite IEAS hypothesis with edge mutliplicity sequence estimated from data
#'   if 'ISA' and deg.hyp = 0: composite ISA hypothesis with edge mutliplicity sequence estimated from data
#' @return
#'  \item{probS}{Probability distribution of Pearson statistic S}
#'  \item{probA}{Probability distribution of divergence statistic A}
#'  \item{summary}{Expected value and variances of test statistics, critical values and significance level
#'  according to asymptotic chi2-distribution, power of tests}
#'  \item{adjusted.stats}{Expected value and variances for adjusted test statistics, preferred adjusted statistics}
#'  \item{adjusted.chi2}{Degrees of freedom for adjusted chi2-distribution}
#'  \item{power.apx}{Power approximations according to adjusted statistics}
#'  \item{degrees.of.freedom}{Degrees of freedom for test performed}
#' @details  to be completed
#' @author Termeh Shafie
#' @seealso [gof_stats],[edge_assignment_probs],[nsumk],[rsm_model]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#' @examples
#' @export
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
