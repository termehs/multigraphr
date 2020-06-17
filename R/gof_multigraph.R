#' @title Goodness of Fit Tests for Random Multigraph Models
#' @description  Goodness of fits tests of multigraph models using Pearson (S) and
#' information divergence (A) test statistics under the random stub matching (RSM)
#' and by independent edge assignments (IEA) model,
#' where the latter is either independent edge assignments of stubs (IEAS) or
#' independent stub assignment (ISA).
#' @param m integer indicating number of edges
#' @param model character string representing assumed model, either RSM, IEAS or ISA
#' @param deg.mod vector of integers with sum 2m representing
#' the modelled degree sequence of the multigraph
#' @param hyp  character string representing testing hypothesis, either IEAS or ISA
#' @param deg.hyp vector of integers with sum 2m representing the hypothetical
#' degree sequence of the multigraph: \cr
#'   - if 'IEAS': simple IEAS hypothesis with fully specified degree sequence deg.hyp\cr
#'   - if 'ISA': simple ISA hypothesis with with fully specified stub assignment probabilities deg.hyp/2m\cr
#'   - if 'IEAS': and deg.hyp = 0: composite IEAS hypothesis with edge mutliplicity sequence estimated from data\cr
#'   - if 'ISA' and deg.hyp = 0: composite ISA hypothesis with edge mutliplicity sequence estimated from data\cr
#' @return
#'  \item{probS}{Probability distribution of Pearson statistic S}
#'  \item{probA}{Probability distribution of divergence statistic A}
#'  \item{summary}{Expected value and variances of test statistics, critical values and significance level
#'  according to asymptotic chi2-distribution, power of tests}
#'  \item{adjusted.stats}{Expected value and variances for adjusted test statistics, preferred adjusted statistics}
#'  \item{adjusted.chi2}{Degrees of freedom for adjusted chi2-distribution}
#'  \item{power.apx}{Power approximations according to adjusted statistics}
#'  \item{degrees.of.freedom}{Degrees of freedom for test performed}
#' @details The tests are performed using goodness-of-fit measures between the
#' edge multiplicity sequence of an observed multigraph,
#' and the expected multiplicity sequence according to a simple or composite hypothesis.
#' Test statistics of Pearson type (S) and of information divergence (A) type are used.
#' @author Termeh Shafie
#' @seealso [gof_stats],[edge_assignment_probs],[nsumk],[rsm_model]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' ## Testing a simple IEAS hypothesis
#' # with degree sequence [6,6,6,2] against
#' # an IEAS model with degree sequence [14,2,2,2]
#' # on a multigrpah with n = 4 nodes and m = 10 edges.
#' deg.mod <- c(14,2,2,2)
#' deg.hyp <- c(6,6,6,2)
#' test1 <- gof_multigraph <- (10, 'IEAS', deg.mod, 'IEAS', deg.hyp)
#'
#' ## Non-null distributiona of test statistics S and A are given by
#' test1$probS
#' test1$probA
#'
#' ## Testing a composite IEAS hypothesis
#' # with degree sequence [15,15,15,15] against
#' # an RSM model with degree sequence [15,15,15,15]
#' # on a multigrpah with n = 4 nodes and m = 30 edges.
#' deg.mod <- c(15,15,15,15)
#' deg.hyp <- c(15,15,15,15)
#' test2 <- gof_multigraph <- (30, 'RSM', deg.mod, 'IEAS', deg.hyp)
#'
#' ## Non-null distributions of test statistics S and A are given by
#' test2$probS
#' test2$probA
#'
#' ## Summary of above tests
#' test1$summary
#' test2$summary
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
