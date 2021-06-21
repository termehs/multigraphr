#'@title Goodness of fit tests for multigraph representations of observed networks
#' @description  Goodness of fits between observed multigraph and specifiedRSM or IEAS hypotheses
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
n <- nrow(adj)

# Observed values = O, Expected values = E
if (type == 'multigraph') {
O <- adj[lower.tri(t(adj), TRUE)]
} else if (type == 'graph') {
adj2 <- adj + diag(adj)
O <- adj2[lower.tri(t(adj2), TRUE)]
}

# expected values depending on whether simple or composite hypothesis
if (deg.hyp == 0) {
  deg.est <- matrix(0,     choose(m+r-1,m), n)
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
      edge_assignment_probs(m, deg.seq, model = hyp)
  }
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
  S <- round(S,3)
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
    # probability distribution of the D values
  D <- round(D, 5) # if you wish to round
  # probability distribution of the A values
  A <- (2 * m * D) / log2(exp(1)) # the asymptotic A statistics
  A <- round(A,3) # if you wish to round
  } else{
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