#'@title Goodness of fit tests for multigraph representations of observed networks
#' @description  Goodness of fits between observed multigraph and specified RSM or IEAS hypotheses
#' using using Pearson (S) and information divergence (A) tests statistics
#' @param m integer giving number of edges in multigraph
#' @param dof  integer giving degrees of freedom of test,
#' r-1 for simple hypotheses and r-n for composite hypotheses where $r = \binom{n+1}{2}$
#' (consider changing these according to preferred adjusted chi-square distributions)
#' @return
#'  \item{summary}{Table including observed test statistics S and A, degrees of freedom for
#'  asymptotic chi-square distribution, and p-value for each test statistics.
#' @details The tests are performed using goodness-of-fit measures between the
#' edge multiplicity sequence of an observed multigraph,
#' and the expected multiplicity sequence according to a simple or composite IEAS hypothesis.
#' Test statistics of Pearson type (S) and of information divergence (A) type are used and summary
#' of test using these two statistics are  given as output.
#' @author Termeh Shafie
#' @seealso [get_edgemultipl_seq],[get_degree_seq],[edge_assignment_probs]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. Journal of Social Structure, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. The Journal of Mathematical Sociology, 40(4), 239-264.
#' @examples
#' ## adjacency matrix of observed network (multigraph), n = 4 nodes , m = 15 edges
#' m_F_PT <- t(matrix(c( 0, 1, 0, 3,
#'                       0, 0, 1, 6,
#'                       0, 0, 0, 3,
#'                       0, 0, 0, 1), nrow= 4, ncol=4))
#' # Testing a simple IEAS hypothesis with degree sequence [4 4 8 14]
#' gof_test(adj, type, 'IEAS', c(4,4,8,14), 9)
#' # Testing a composite IEAS hypothesis
#' gof_test(adj, type, 'IEAS', 0, 6)
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
}

# expected values depending on whether simple or composite hypothesis
if (sum(deg.hyp) == 0) {
  deg.est <- get_degree_seq(adj, type)
  Q.seq <-  edge_assignment_probs(m, deg.est, hyp)
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
Q.seq <- edge_assignment_probs(m, deg.hyp, hyp)
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
