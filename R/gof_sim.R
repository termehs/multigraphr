#' @title Goodness of fit test simulations for random multigraph models
#' @description  Goodness of fits test simulations of specified multigraph models using Pearson (\emph{S}) and
#' information divergence (\emph{A}) test statistics under the random stub matching (RSM)
#' and the independent edge assignments (IEA) model,
#' where the latter is either independent edge assignments of stubs (IEAS) or
#' independent stub assignment (ISA).
#'
#' These can be used to check the reliability of the tests by examining the exact probability
#' distributions of the test statistics and their fit to their asymptotic  chi^2-distributions.
#' Only practical for small multigraphs as exact distributions are calculated.
#' @importFrom stats dmultinom
#' @param m integer giving number of edges in multigraph.
#' @param model character string representing assumed model, either \code{'RSM'}, \code{'IEAS'} or \code{'ISA'}.
#' @param deg.mod vector of integers with the sum equal to 2\code{m} representing
#' the degree sequence of the multigraph under specified model.
#' @param hyp  character string representing the hypothesized null model, either \code{'IEAS'} or \code{'ISA'}.
#' @param deg.hyp vector of integers with the sum equal to to 2\code{m} representing the hypothetical
#' degree sequence of the multigraph under the null model: \cr
#'   - if \code{hyp = 'IEAS'}, then simple IEAS hypothesis with fully specified degree sequence \code{deg.hyp}\cr
#'   - if \code{hyp = 'ISA'}, then simple ISA hypothesis with with fully specified stub assignment probabilities \code{deg.hyp}/2\code{m}\cr
#'   - if \code{hyp = 'IEAS'} and \code{deg.hyp = 0}, then composite IEAS hypothesis with edge multiplicity sequence estimated from data\cr
#'   - if \code{hyp = 'ISA'} and \code{deg.hyp = 0}, then composite ISA hypothesis with edge multiplicity sequence estimated from data\cr
#' @return Output is generated from function \code{\link{gof_stats}}:
#'  \item{test.summary}{Expected value and variances of test statistics (\code{stat}),
#'  critical values (\code{cv}) according to asymptotic chi^2-distribution and
#'  according to cdf's of test statistics,
#'  significance level (alpha) according to asymptotic chi^2 distribution,
#'  power of tests (\code{P(stat>cv)}), critical values and power
#'  according to the distributions of test statistics (\code{cv(stat)}
#'  and \code{ P(Stat>cv(Stat))}).}
#'  \item{degrees.of.freedom}{Degrees of freedom for tests performed.}
#'  \item{probS}{Probability distributions of Pearson statistic \code{S}.}
#'  \item{probA}{Probability distributions of information divergence statistic \code{A}.}
#'  \item{adjusted.stats}{Expected values and variances for adjusted test statistics,
#'  preferred adjusted statistics.}
#'  \item{adjusted.chi2}{Degrees of freedom for adjusted  chi^2-distribution.}
#'  \item{power.apx}{Power approximations according to adjusted statistics.}
#' @details The tests are performed using goodness-of-fit measures between simulated
#' edge multiplicity sequence of a multigraph according to an RSM or IEA model,
#' and the expected multiplicity sequence according to a simple or composite IEA hypothesis.
#'
#' Test statistics of Pearson type (\emph{S}) and
#' of information divergence (\emph{A}) type are used and summary
#' of tests given these two statistics are given as output. The adjusted statistics and
#' chi^2-distributions are useful for better power calculations.
#' @author Termeh Shafie
#' @seealso \code{\link{gof_stats}},\code{\link{get_edge_assignment_probs}},
#' \code{\link{nsumk}},\code{\link{rsm_model}}
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' # Testing a simple IEAS hypothesis with degree sequence [6,6,6] against
#' # an IEAS model with degree sequence [14,2,2] on a multigraph with n = 3 nodes and m = 9 edges
#' deg.mod <- c(14,2,2)
#' deg.hyp <- c(6,6,6)
#' test1 <- gof_sim(9, 'IEAS', deg.mod, 'IEAS', deg.hyp)
#'
#' # Non-null distributions (pdf's and cdf's) of test statistics S and A are given by
#' test1$probS
#' test1$probA
#'
#' # Testing a simple ISA hypothesis with degree sequence [6,6,6] against
#' # an IEAS model with degree sequence [12,3,3] on a multigraph with n = 3 nodes and m = 9 edges
#' deg.mod <- c(12,3,3)
#' deg.hyp <- c(6,6,6)
#' test2 <- gof_sim(9, 'IEAS', deg.mod, 'ISA', deg.hyp)
#'
#' # Testing a composite IEAS hypothesis against
#' # an RSM model with degree sequence [6,6,6,6] on a multigraph with n = 4 nodes and m = 20 edges.
#' deg.mod <- c(6,6,6,6)
#' test3 <- gof_sim(12, 'RSM', deg.mod, 'IEAS', deg.hyp = 0)
#'
#' @export
#'
gof_sim <- function(m, model, deg.mod, hyp, deg.hyp) {
  n <- length(deg.mod)
  r <- choose(n + 1, 2)

  if (sum(deg.mod)/2 != m)
    stop("number of edges must be half the sum of the degree sequence")

  # model specification: IEAS or ISA
  if (model == 'IEAS') {
    Q.seq <- get_edge_assignment_probs(m, deg.mod, model = 'IEAS')
    # use IEA model to find probabilities of each multiplicity sequence/multigraph
    m.seq <- nsumk(r, m)
    prob.mg <-
      apply(m.seq, 1, function(x)
        dmultinom(x, prob = Q.seq))
  }
  else if (model == 'ISA') {
    Q.seq <- get_edge_assignment_probs(m, deg.mod, model = 'ISA')
    # use ISA model to find probabilities of each multiplicity sequence/multigraph
    m.seq <- nsumk(r, m)
    prob.mg <-
      apply(m.seq, 1, function(x)
        dmultinom(x, prob = Q.seq))
  }
  else if (model == 'RSM') {
    # use RSM model to find probabilities of each multiplicity sequence/multigraph
    rsm <- rsm_model(deg.mod)
    m.seq <- as.matrix(rsm$m.seq)
    prob.mg <- rsm$prob.dists[, 1]
  }

  # according to hypothesis, get the values of test statistics
  # hypothesis IEAS
  if (hyp == 'IEAS') {
    if (sum(deg.hyp) > 0) {
      if (sum(deg.hyp)/2 != m)
        stop("number of edges must be half the sum of the degree sequence")
      dof <- r - 1
      # fully specified with deg.hyp
      Q.seq <- get_edge_assignment_probs(m, deg.hyp, model = 'IEAS')
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
    else if (sum(deg.hyp) == 0) {
      dof <- r - n
      # estimated deg.est from data
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
          get_edge_assignment_probs(m, deg.seq, model = 'IEAS')
      }
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
  }

  # hypothesis ISA
  else if (hyp == 'ISA') {
    if (sum(deg.hyp) > 0) {
      if (sum(deg.hyp)/2 != m)
        stop("number of edges must be half the sum of the degree sequence")
      dof <- r - 1
      # fully specified with deg.hyp
      Q.seq <- get_edge_assignment_probs(m, deg.hyp, model = 'ISA')
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
    else if (sum(deg.hyp) == 0) {
      dof <- r - n
      # estimated deg.est from data
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
          get_edge_assignment_probs(m, deg.seq, model = 'ISA')
      }
      moms <- gof_stats(m, dof, m.seq, prob.mg, Q.seq)
    }
  }
}
