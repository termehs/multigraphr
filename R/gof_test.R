#'@title Goodness of fit tests
#' @description  Goodness of fit tests between observed edge multiplicity sequence and expected edge multiplicity sequence given a
#' specified RSM or IEAS hypotheses using using Pearson (S) and information divergence (A) tests statistics
#' @param adj Matrix of integer representing graph adjacency matrix
#' @param type Equals 'graph' if adjacency matrix is for graphs (default)
#' @param hyp  character string representing the hypothesized model (null), either IEAS or ISA
#' @param deg.hyp vector of integers with sum equal to 2m representing the null
#' degree sequence of the multigraph: \cr
#'   - if 'IEAS': simple IEAS hypothesis with fully specified degree sequence deg.hyp\cr
#'   - if 'ISA': simple ISA hypothesis with with fully specified stub assignment probabilities deg.hyp/2m\cr
#'   - if 'IEAS': and deg.hyp = 0: composite IEAS hypothesis with edge multiplicity sequence estimated from data\cr
#'   - if 'ISA' and deg.hyp = 0: composite ISA hypothesis with edge multiplicity sequence estimated from data\cr
#' @param m integer giving number of edges in multigraph
#' @param dof  integer giving degrees of freedom of test,
#' r-1 for simple hypotheses and r-n for composite hypotheses where $r = \binom{n+1}{2}$
#' @return
#'  \item{summary} Table including observed test statistics S and A, degrees of freedom for
#'  their asymptotic  χ²-distribution, and p-values for tests performed
#' @details This function can be used to test whether there is a significant difference between
#'  observed multigraph and the expected multiplicity sequence according
#'  to a simple or composite IEAS hypothesis.
#'  Test statistics of Pearson (S) and of information divergence (A) type are used and
#'  test summary based on these two statistics are given as output.
#'  \emph{P}-values indicate whether the null we have sufficient evidence to reject the null
#'  that there is no significant difference between the observed and expected edge multiplicity sequence.
#' @author Termeh Shafie
#' @seealso [gof_sim],[get_degree_seq],[nsumk],[rsm_model],[gof_stats]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. Journal of Social Structure, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. The Journal of Mathematical Sociology, 40(4), 239-264.
#' @examples
#' ## adjacency matrix of observed network (multigraph), n = 4 nodes , m = 15 edges
#' adj <- t(matrix(c( 0, 1, 0, 3,
#'                    0, 0, 1, 7,
#'                    0, 1, 0, 3,
#'                    3, 6, 3, 2), nrow= 4, ncol=4))
#' deg <- get_degree_seq(adj, 'multigraph')
#' # Testing a simple IEAS hypothesis with above degree sequence
#' gof_test(adj, type = 'multigraph', 'IEAS', deg, 9)
#' # Testing a composite IEAS hypothesis
#' gof_test(adj, type  = 'multigraph', 'IEAS', 0, 6)
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
