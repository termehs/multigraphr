#' @title Exact probability distributions and moments of goodness of fit statistics
#' @description  Goodness of fit between two specified edge multiplicity sequences
#' (e.g. observed vs. expected).
#' Pearson (\emph{S}) and information divergence (\emph{A}) tests statistics are used and
#' the exact distribution of these statistics,  their asymptotic chi^2-distributions,
#' and their first two central moments are calculated using this function.
#' Only practical for small multigraphs.
#' @importFrom stats pchisq
#' @param m integer giving number of edges in multigraph.
#' @param dof  integer giving degrees of freedom of test performed.
#' @param m.seq  matrix of integers, each row representing the
#' edge multiplicity sequence of a multigraph (which correspond to observed values).
#' @param prob.mg  numerical vector representing a given probability distribution of
#' multigraphs/edge multiplicity sequences in \code{m.seq}.
#' @param Q.seq  a numeric vector representing the hypothetical edge assignment probabilities
#' to all possible vertex pair sites (from which expected values are calculate).
#' @return
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
#' @details The tests are performed using goodness-of-fit measures between two edge multiplicity sequences
#' (e.g. observed vs. expected).
#'
#' Test statistics of Pearson type (\emph{S}) and
#' of information divergence (\emph{A}) type are used and summary
#' of tests given these two statistics are given as output. The adjusted statistics and
#' chi^2-distributions are useful for better power calculations.
#' @author Termeh Shafie
#' @seealso \code{\link{gof_sim}},\code{\link{get_edge_assignment_probs}},
#' \code{\link{nsumk}}
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' # Generate a set of edge multiplicity sequences (random multigraphs) and
#' # its probability distribution using rsm_model() with degree sequence [4,4,6,6]
#' rsm <- rsm_model(deg.seq = c(4,4,6,6))
#' mg <- as.matrix(rsm$m.seq)
#' mg.p <- rsm$prob.dists[, 1]
#'
#' # Generate edge assignment probabilities from which the second set of
#' # edge multiplicity sequences is generated from using the iea_model()
#' deg.f <- (4*5)/2 - 1
#' eap <- get_edge_assignment_probs(m = 10,
#'                    deg.seq = c(4,4,6,6), model = 'IEAS')
#'
#' # Perform the test
#' test <- gof_stats(m = 10, dof = deg.f,
#'                    m.seq = mg, prob.mg = mg.p, eap)
#'
#' @export
#'
gof_stats <- function(m, dof, m.seq, prob.mg, Q.seq) {
  # the chi-square distribution with dof degrees of freedom, cv = critical value
  cv <- dof + sqrt(8 * dof)
  prob.Chi <- pchisq(cv, dof) # prob(chi(dof) < cv)
  alpha <- 1 - prob.Chi # prob(chi(dof) > cv )


  # Observed values = O, Expected values = E
  O = m.seq
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
  S <-
    round(S, 5) # if you wish to round before finding unique values to make outcome space smaller
  S.uni <- sort(unique(S))
  # probability distribution of the S values
  prob.tmp <- vector()
  prob.S <- vector()
  for (i in 1:length(S.uni)) {
    for (j in 1:length(S)) {
      if (S.uni[i] == S[j]) {
        prob.tmp <- c(prob.tmp, prob.mg[j])
      }
    }
    prob.S[i] <- sum(prob.tmp)
    prob.tmp <- vector()
  }
  cumprob.S <- cumsum(prob.S)
  ExpS <- sum(S.uni * prob.S)
  VarS <- sum((S.uni ^ 2 * prob.S)) - ExpS ^ 2
  # prob(Shat > cv)
  S.g.cv <- sum(prob.S[S.uni > cv])
  # prob(Shat > cv)
  cvS <- ExpS + 2 * sqrt(VarS)
  S.g.cvS <- sum(prob.S[S.uni > cvS])
  S.out <-
    round(cbind(ExpS, VarS, cv, alpha, S.g.cv, cvS, S.g.cvS), 5) # for output
  # approximate statistics (adjusted chi square): S' and S''
  S.prim <- S.uni * (floor(ExpS) / ExpS)
  ExpS.prim <- ExpS
  VarS.prim <- (2 * ExpS ^ 2) / (floor(ExpS))
  S.bis <- S.uni * (dof / ExpS)
  ExpS.bis <- ExpS
  VarS.bis <- (2 * ExpS ^ 2) / dof

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
  # A <- round(A,3) # if you wish to round
  A <- (2 * m * D) / log2(exp(1)) # the asymptotic A statistics
  A.uni <- sort(unique(A))
  prob.tmp <- vector()
  prob.A <- vector()
  for (i in 1:length(A.uni)) {
    for (j in 1:length(A)) {
      if (A.uni[i] == A[j]) {
        prob.tmp <- c(prob.tmp, prob.mg[j])
      }
    }
    prob.A[i] <- sum(prob.tmp)
    prob.tmp <- vector()
  }
  cumprob.A <- cumsum(prob.A)
  ExpA <- sum(A.uni * prob.A)
  VarA <- sum((A.uni ^ 2 * prob.A)) - ExpA ^ 2
  # prob(Shat > cv)
  A.g.cv <- sum(prob.A[A.uni > cv])

  # prob(Shat > cv)
  cvA <- ExpA + 2 * sqrt(VarA)
  A.g.cvA <- sum(prob.A[A.uni > cvA])
  A.out <-
    round(cbind(ExpA, VarA, cv, alpha, A.g.cv, cvA, A.g.cvA), 5)

  # approximate statistics (adjusted chi square): S' and S''
  A.prim <- A.uni * (floor(ExpA) / ExpA)
  ExpA.prim <- ExpA
  VarA.prim <- (2 * ExpA ^ 2) / (floor(ExpA))
  A.bis <- A.uni * (dof / ExpA)
  ExpA.bis <- ExpA
  VarA.bis <- (2 * ExpA ^ 2) / dof

  # output 1: probability distributions of S and A
  prob.Sout <-
    as.data.frame(round(cbind(S.uni, prob.S, cumprob.S), 5))
  colnames(prob.Sout) <- c('S=s', 'P(S=s)', 'P(S<s)')
  prob.Aout <-
    as.data.frame(round(cbind(A.uni, prob.A, cumprob.A), 5))
  colnames(prob.Aout) <- c('A=a', 'P(A=a)', 'P(A<a)')
  # output 2: summary for statistics S and A
  gof.sum <- as.data.frame(rbind(S.out, A.out))
  stat <- c('S', 'A')
  gof.sum <- cbind(stat, gof.sum)
  colnames(gof.sum) <-
    c('Stat',
      'E(Stat)',
      'V(Stat)',
      'cv',
      'alpha',
      'P(Stat>cv)',
      'cv(Stat)',
      'P(Stat>cv(Stat))')

  # output 3: moments of approximate statistics (just difference in variance)
  Exp.out <-
    cbind(ExpS, ExpS.prim, ExpS.bis, ExpA, ExpA.prim, ExpA.bis)
  Var.out <-
    cbind(VarS, VarS.prim, VarS.bis, VarA, VarA.prim, VarA.bis)

  # Determine the best S and A approximations
  Best.apx <- rep(0, 6)
  var.adj.S <- Var.out[1, 2:3] # variance of S' and S''
  var.S <- Var.out[1, 1] # variance of S
  if (Var.out[2] != Var.out[3]) {
    idx <- which.min(abs(var.adj.S - var.S))
    if (idx == 1) {
      Best.apx[idx + 1] <- 1
    } else if (idx == 2) {
      Best.apx[idx + 1] <- 1
    }
  }
  var.adj.A <- Var.out[1, 5:6] # variance of A' and A''
  var.A <- Var.out[1, 4] # variance of A
  if (Var.out[5] != Var.out[6]) {
    idx <- which.min(abs(var.adj.A - var.A))
    if (idx == 1) {
      Best.apx[idx + 4] <- 1
    } else if (idx == 2) {
      Best.apx[idx + 4] <- 1
    }
  }
  adj.out <- as.data.frame(rbind(Exp.out, Var.out, Best.apx))
  row.names(adj.out) <- c('Exp', 'Var', 'Best Adj.')
  colnames(adj.out) <-
    c('S', 'Sprim', 'Sbis', 'A', 'Aprim', 'Abis')

  # output 4: adjusted chi square distributions (their degrees of freedom)
  S.comp <- 2*(ExpS^2)/var.S
  if (abs(S.comp - dof) < abs(S.comp - floor(ExpS))) {
    adj.chi2.S <- dof
  } else {
    adj.chi2.S <- floor(ExpS)
  }
  A.comp<- 2*(ExpA^2)/var.A
  if (abs(A.comp - dof) < abs(A.comp - floor(ExpA))) {
    adj.chi2.A <- dof
  } else {
    adj.chi2.A <- floor(ExpA)
  }
  adj.chi2 <- as.data.frame(cbind(adj.chi2.S, adj.chi2.A))
  colnames(adj.chi2) <- c('df(S)', 'df(A)')

  # output 5: power approximations with adjusted S and A statistics
  tmp <- cv * floor(ExpS) / ExpS
  powerS.prim <- pchisq(tmp, floor(ExpS))
  tmp2 <- cv * (2 * ExpS ^ 2 / dof) / ExpS
  powerS.bis <- pchisq(tmp2, floor(ExpS))
  tmp <- cv * floor(ExpA) / ExpA
  powerA.prim <- pchisq(tmp, floor(ExpA))
  tmp2 <- cv * (2 * ExpA ^ 2 / dof) / ExpA
  powerA.bis <- pchisq(tmp2, floor(ExpA))
  power.apx <- as.data.frame(round(cbind(
    powerS.prim, powerS.bis, powerA.prim, powerA.bis
  ), 5))
  colnames(power.apx) <- c('Sprim', 'Sbis', 'Aprim', 'Abis')

  # output list
  listout <-
    list(
      "test.summmary" = gof.sum,
      "degrees.of.freedom" = dof,
      "probS" = prob.Sout,
      "probA" = prob.Aout,
      "adjusted.stats" = adj.out,
      "adjusted.chi2" = adj.chi2,
      "power.apx" = power.apx
    )
  return(listout)
}
