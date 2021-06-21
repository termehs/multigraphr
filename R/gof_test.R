#'@title Goodness of fit tests for multigraph representations of observed networks
#' @description  Goodness of fits between observed multigraph and specifiedRSM or IEAS hypotheses
#' using using Pearson (S) and information divergence (A) tests statistics
#' @param m integer giving number of edges in multigraph
#' @param dof  integer giving degrees of freedom of test (consider changing these accpording to 
#' preferred adjuested chi-square distributions)
#' @param m.seq  vector of integers representing observed multigraph
#' @param Q.seq  A numeric vector representing the edge assignment probabilities
#' to all possible vertex pair sites according to a simple or composite IEAS hypothesis
#' @return
#'  \item{test.summary}{Expected value and variances of test statistics (stat),
#'  critical values (cv) according to asymptotic chi2 distribution and
#'  according to cdf's of test statistics,
#'  significance level (alpha) according to
#'  asymptotic chi2 distribution, and power of tests (P(stat>cv))}
#' @details The tests are performed using goodness-of-fit measures between the
#' edge multiplicity sequence of an observed multigraph,
#' and the expected multiplicity sequence according to a simple or composite IEAS hypothesis.
#' Test statistics of Pearson type (S) and of information divergence (A) type are used and summary
#' of tests given these two statistics are given as output.
#' @author Termeh Shafie
#' @seealso [get_edgemultipl_seq],[edge_assignment_probs]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. Journal of Social Structure, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. The Journal of Mathematical Sociology, 40(4), 239-264.
#' @examples
#' ## For examples see gof_multigraph()
#' @export