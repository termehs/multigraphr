#'@title Goodness of fit tests for random multigraph models
#' @description  Goodness of fit tests between an observed edge multiplicity sequence and
#' an expected edge multiplicity sequence according to specified RSM or IEA hypotheses
#' using Pearson (\emph{S}) and information divergence (\emph{A}) tests statistics.
#' @importFrom stats pchisq
#' @param adj matrix of integer representing graph adjacency matrix.
#' @param type equals \code{'graph'} if adjacency matrix is for graphs (default),
#' equals \code{'multigraph'} if it is the equivalence of the adjacency matrix for multigraphs
#' (with matrix diagonal representing loops double counted).
#' @param hyp  character string representing the null model, either \code{'IEAS'} or \code{'ISA'}.
#' @param deg.hyp vector of integers with the sum equal to to 2\code{m} representing the
#' degree sequence of the multigraph under the null model: \cr
#'   - if \code{hyp = 'IEAS'}, then simple IEAS hypothesis with fully specified degree sequence \code{deg.hyp}\cr
#'   - if \code{hyp = 'ISA'}, then simple ISA hypothesis with with fully specified stub assignment probabilities \code{deg.hyp}/2\code{m}\cr
#'   - if \code{hyp = 'IEAS'} and \code{deg.hyp = 0}, then composite IEAS hypothesis with edge multiplicity sequence estimated from data\cr
#'   - if \code{hyp = 'ISA'} and \code{deg.hyp = 0}, then composite ISA hypothesis with edge multiplicity sequence estimated from data\cr
#' @param dof  integer giving degrees of freedom of test,
#' \emph{r-1} for simple hypotheses and \emph{r-n} for composite hypotheses where \emph{r = n(n+1)/2}
#' @return
#'  \item{summary}{Data frame including observed  values of test statistics \code{S} and \code{A},
#'  degrees of freedom for the asymptotic chi^2-distributions of tests statistics,
#'  and \emph{p}-values for the tests performed.}
#' @details This function can be used to test whether there is a significant difference between
#'  observed multigraph and the expected multiplicity sequence according
#'  to a simple or composite IEA hypothesis.
#'
#'  Test statistics of Pearson (\emph{S}) and of information divergence (\emph{A}) type are used and
#'  test summary based on these two statistics are given as output.
#'
#'  \emph{p}-values indicate whether we have sufficient evidence to reject the null
#'  that there is no significant difference between the observed and
#'  expected edge multiplicity sequence.
#' @author Termeh Shafie
#' @seealso \code{\link{get_degree_seq}},\code{\link{get_edge_assignment_probs}},
#' \code{\link{gof_sim}} to check the reliability of your test
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' # Adjacency matrix of observed network (multigraph), n = 4 nodes , m = 15 edges
#' A <- t(matrix(c( 0, 1, 0, 3,
#'                    0, 0, 1, 7,
#'                    0, 1, 0, 3,
#'                    3, 6, 3, 2), nrow= 4, ncol=4))
#' deg <- get_degree_seq(adj = A, type = 'multigraph')
#'
#' # Testing a simple IEAS hypothesis with above degree sequence
#' gof_test(adj = A, type = 'multigraph', hyp = 'IEAS', deg.hyp = deg, dof = 9)
#'
#' # Testing a composite IEAS hypothesis
#' gof_test(adj = A, type  = 'multigraph', hyp = 'IEAS', deg.hyp = 0, dof = 6)
#' @export
#'
gof_test <- function(adj, type, hyp, deg.hyp, dof) {
m <- sum(adj) / 2
n <- nrow(adj)

# Observed values = O, Expected values = E
if (type == 'multigraph') {
O <- adj[lower.tri(t(adj), TRUE)]
} else if (type == 'graph') {
adj2 <- adj + diag(adj)
O <- adj2[lower.tri(t(adj2), TRUE)]
} else {
  stop("type must be defined as either 'graph' or 'multigraph'")
}

# expected values depending on whether simple or composite hypothesis
if (sum(deg.hyp) == 0) {
  deg.est <- get_degree_seq(adj, type)
  Q.seq <-  get_edge_assignment_probs(m, deg.est, hyp)
  E = m * Q.seq
  # S = Pearson goodness-of-fit statistic
  S = (t(O) - E) ^ 2 / E
  S[is.infinite(S) | is.na(S)] = 0
  S = round(sum(S),3)
  # D = divergence goodness-of-fit statistic
  D <- (O / m) * log2(O / E)
  D[is.infinite(D) | is.na(D)] = 0
  A <- (2 * m * D) / log2(exp(1))
  A <- round(sum(A),3)
  } else if (sum(deg.hyp) > 0) {
  Q.seq <- get_edge_assignment_probs(m, deg.hyp, hyp)
  E = m * Q.seq
# S = Pearson goodness-of-fit statistic
  S = (t(O) - E) ^ 2 / E
  S[is.infinite(S) | is.na(S)] = 0
  S = round(sum(S),3)
# D = divergence goodness-of-fit statistic
  D <- (O / m) * log2(O / E)
  D[is.infinite(D) | is.na(D)] = 0
  A <- (2 * m * D) / log2(exp(1))
  A <- round(sum(A),3)
}
# p-values
pval_S <- round(1 - pchisq(S, dof),3) # prob(chi(dof) < S)
pval_A <- round(1 - pchisq(A, dof),3) # prob(chi(dof) < A)

# output: test summary using S and A
S.out <- cbind(S, pval_S)
A.out <- cbind(A, pval_A)
test.sum <- rbind(S.out, A.out)
stat <- c('S', 'A')
summary <- as.data.frame(cbind(stat, dof, test.sum))
colnames(summary) <- c('Stat', 'dof', 'Stat(obs)',  'p-value')

return(summary)

 }
