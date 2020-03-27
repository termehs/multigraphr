#' @title Moments of Goodness of Fit Statistics
#' @description
#' @param
#' @return
#' @details  To be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#' @examples
#'
gof_stats <- function(m, dof, m.seq, prob.mg, Q.seq) {
  # the chi-square distribution with dof degrees of freedom, cv = critical value
  cv <- dof + sqrt(8 * dof)
  prob.Chi <- pchisq(cv, dof) # prob(chi(dof) < cv)
  alpha <- 1 - prob.Chi # prob(chi(dof) > cv )
  # S = Pearson goodness-of-fit statistic
  # Observed values = O, Expected values = E
  O = m.seq
  E = m * Q.seq
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
  #] D.uni <- sort(unique(D))
  # prob.tmp <- vector()
  # prob.D <- vector()
  # for(i in 1:length(D.uni)){
  #   for(j in 1:length(D)){
  #     if(D.uni[i] == D[j]){
  #       prob.tmp <- c(prob.tmp, prob.mg[j])
  #     }
  #   }
  #   prob.D[i] <- sum(prob.tmp)
  #   prob.tmp <- vector()
  # }
  # cumprob.D <- cumsum(prob.D)
  # ExpD <- sum(D.uni*prob.D)
  # VarD <- sum((D.uni^2*prob.D))-ExpD^2
  # # prob(Shat > cv)
  # D.g.cv <- sum(prob.D[D.uni > cv])
  # # prob(Shat > cv)
  # cvD <- ExpD+2*sqrt(VarD)

  # probability distribution of the A values
  # A <- round(A,3) # if you wish to round
  A <- (2 * m * D) / log2(exp(1)) # the asymptotic T statistics
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
    c('stat',
      'exp',
      'var',
      'cv',
      'alpha',
      'stat>cv',
      'cv(stat)',
      'stat>cv(stat)')


  # output 3: moments of approximate statistics (just differene in variance)
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

  # output 4: power approximations with S and A
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

  # output 5: plot the distributions of S and A, together with Chi
  listout <-
    list(
      "probS" = prob.Sout,
      "probA" = prob.Aout,
      "summmary" = gof.sum,
      "adjusted.stats" = adj.out,
      "power.apx" = power.apx,
      "degrees.of.freedom" = dof
    )
  return(listout)
}
