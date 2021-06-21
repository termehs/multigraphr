#'@title Goodness of fit tests for multigraph representations of observed networks
#' @description  Goodness of fits between observed multigraph and specifiedRSM or IEAS hypotheses
#' using using Pearson (S) and information divergence (A) tests statistics
#' @param m integer giving number of edges in multigraph
#' @param dof  integer giving degrees of freedom of test, 
#' r-1 for simple hypotheses and r-n for composite hypotheses where $r = \binom{n+1}{2}$
#' (consider changing these according to preferred adjusted chi-square distributions)
#' @param m.seq  vector of integers representing observed multigraph
#' @param Q.seq  A numeric vector representing the edge assignment probabilities
#' to all possible vertex pair sites according to a simple or composite IEAS hypothesis
#' @return
#'  \item{summary}{Table including observed test statistics S and A, degrees of freedom for
#'  asymptotic chi-square distribution, and p-value for each test statistics. 
#' @details The tests are performed using goodness-of-fit measures between the
#' edge multiplicity sequence of an observed multigraph,
#' and the expected multiplicity sequence according to a simple or composite IEAS hypothesis.
#' Test statistics of Pearson type (S) and of information divergence (A) type are used and summary
#' of test using these two statistics are  given as output.
#' @author Termeh Shafie
#' @seealso [get_edgemultipl_seq],[edge_assignment_probs]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. Journal of Social Structure, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. The Journal of Mathematical Sociology, 40(4), 239-264.
#' @examples
#' ## For examples see gof_multigraph() 
#' @export
#' 
gof_test <- function(adj, type, hyp, deg.hyp, dof) {
m <- sum(adj) / 2

# Observed values = O, Expected values = E
O <- get_edgemultip_seq(adj, type)  
Q.seq <- edge_assignment_probs(m, deg.hyp, hyp)
E = m * Q.seq
# S = Pearson goodness-of-fit statistic
if (is.vector(E)) {
  S = (t(O) - E) ^ 2 / E
  S[is.infinite(S) | is.na(S)] = 0
  S = colSums(S)
} else{
  S = (O - E) ^ 2 / E
  S[is.infinite(S) | is.na(S)] = 0
  S = rowSums(S)
}
# D = divergence goodness-of-fit statistic
if (is.vector(E)) {
  D <- (O / m) * log2(O %*% diag(1 / E))
  D[is.infinite(D) | is.na(D)] = 0
  D <- rowSums(D)
} else{
  D <- (O / m) * log2(O * (1 / E))
  D[is.infinite(D) | is.na(D)] = 0
  D <- rowSums(D)
}
# p-values
pval_S <- 1 - pchisq(S, dof) # prob(chi(dof) < S)
pval_S <- 1 - pchisq(S, dof) # prob(chi(dof) < A)

# output: test summary using S and A
S.out <- cbind(S, pval_S)
A.out <- cbind(S, pval_A)
test.sum <- rbind(S.out, A.out)
stat <- c('S', 'A')
summary <- cbind(stat, dof, test.sum)
colnames(summary) <- c('Stat', 'dof', 'p-value = P(Stat < chi2(dof)')

return(summary)
 }