#' @title Random stub matching model for multigraphs
#' @description  Given a specified degree sequence,
#' this function finds all unique multigraphs represented
#' by their edge multiplicity sequences. Different
#' complexity statistics together with their
#' probability distributions and moments are calculated.
#' @param deg.seq vector of integers representing the degree sequence of a multigraph.
#' @return
#' \item{m.seq}{possible multigraphs represented by edge multiplicity sequences}
#' \item{prob.dists}{probability distribution of the multigraphs/edge multiplicity sequences,
#' and the probability distributions of the statistics \emph{number of loops},
#' \emph{number of multiple edges}, and
#' \emph{simple graph} (logical) for each multigraph}
#' \item{M}{summary of moments and interval estimates for
#' \emph{number of loops} and \emph{number of multiple edges} (\code{M1} and \code{M2})).}
#' @details  The probability distributions of all unique multigraphs given fixed degree sequence,
#' together with the first two central moments and interval estimates of the statistics
#' \emph{M1 = number of loops} and \emph{M2 = number of multiple edges}, under the RSM model are calculated.
#'
#' For other structural statistics and for large multigraphs,
#' use the IEA approximation of the RSM model via the function \code{\link{iea_model}}
#' @author Termeh Shafie
#' @seealso \code{\link{get_degree_seq}}, \code{\link{get_edge_multip_seq}}, \code{\link{iea_model}}
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' # Given a specified degree sequence
#' D <- c(2,2,3,3) # degree sequence
#' mod1 <- rsm_model(D)
#' mod1$m.seq
#' mod1$prob.dists
#' mod1$M
#'
#' # Given an observed graph
#' A <-  matrix(c(0, 1, 2,
#'                 1, 2, 1,
#'                 2, 1, 2), nrow=3, ncol=3)
#' D <-  get_degree_seq(adj = A, type = 'graph')
#' mod2 <- rsm_model(D)
#' @export
#'
rsm_model <- function(deg.seq) {
  if (sum(deg.seq) %% 2 == 1)
    stop("sum of degree sequence must be an even number")
  m <- sum(deg.seq / 2)
  n <- length(deg.seq)
  k <- m - 1
  r <- choose(n + 1, 2)

  #initial edge sequence (read as labelled 2-tuples of connected nodes)
  s <- vector()
  edge.seq <- vector()
  for (i in 1:n) {
    s <- rep(i, deg.seq[i])
    edge.seq <- c(edge.seq, s)
  }
  s <- edge.seq

  # inital edge list
  z <- matrix(0, m, 2)
  for (i in 1:m) {
    z[i, 1] = edge.seq[2 * i - 1]
    z[i, 2] = edge.seq[2 * i]
  }

  # create all possible sequences of edges given degree sequence
  while (k > 0)
  {
    W <- vector("integer", 0)
    for (i in k:m) {
      Wz <- z[i,]
      W <- c(W, Wz)
    }
    tmp <- c(W[1], W[2] + 1)
    W <- sort(W)
    if (length(which(tmp[1] == W)) <= length(which(tmp[2] <= W))) {
      i1 <- which(tmp[1] == W)
      i2 <- which(tmp[2] <= W)
      W1 <- vector("integer", 0)
      del <- vector("integer", 0)
      for (i in 1:min(length(i1), length(i2))) {
        w <- c(W[i1[i]], W[i2[i]])
        W1 <- c(W1, w)
        del <- c(del, c(i1[i], i2[i]))
      }
      W <- W[-del]
      ss <- c(W1, W)
      s[(2 * m - length(ss) + 1):(2 * m)] <- ss
      edge.seq <- rbind(edge.seq, s, deparse.level = 0)
      for (i in 1:m) {
        z[i, 1] = s[2 * i - 1]
        z[i, 2] = s[2 * i]
      }
      k = m - 1
    } else {
      k = k - 1
    }
  }

  # for each possible edge sequence/multigraph given degree sequence, count number of loops and number of multiple edges using:
  # edge multiplicity sequence = m.seq
  # ordered edge multiplicity sequence = m.seq.star
  m.seq <- vector()
  m1 <- vector()
  m2 <- vector()
  m3 <- vector()
  m4 <- vector() # m1+m3 which is equal to 0 iff graph is simple
  edge.perm <- vector()
  mz <- vector()
  tz <- vector()
  t <- vector()
  shifts <- vector()
  tot <- vector()
  simple <- vector()
  for (g in 1:nrow(edge.seq)) {
    for (i in 1:m) {
      z[i, 1] = edge.seq[g, 2 * i - 1]
      z[i, 2] = edge.seq[g, 2 * i]
    }
    A <- matrix(0, n, n) #adjacency matrix
    f <- z
    f[, c(1, 2)] <- z[, c(2, 1)]
    for (j in 1:m) {
      b <- vector()
      c <- vector()
      for (i in 1:m) {
        b <- rbind(b, as.integer(z[j, ] == z[i, ])) #for multiple edges
        c <- rbind(c, as.integer(z[j, ] == f[i, ])) #for loops
      }
      if (length(rowSums(b) == 2) > 1) {
        A[z[j, 1], z[j, 2]] <- length(which(rowSums(b) == 2))
      } else {
        A[z[j, 1], z[j, 2]] <- 1
      }
    }

    # number of loops m1
    m1[g] <- sum(diag(A))

    # number of edges that are not loops m2
    m2[g] <- m - m1[g]

    # number of multiple edges m3
    m3[g] <- sum(A[which(A > 1)] - 1)

    # possible permutations in each edge sequence
    if (m3[g] > 0) {
      mult <- A[which(A > 1)]
      edge.perm[g] <- factorial(m) / (prod(factorial(mult)))
    } else {
      edge.perm[g] = factorial(m)
    }

    # the edge multiplicity sequence m.seq
    tmp <- t(A)[lower.tri(t(A), diag = TRUE)]
    m.seq <- rbind(m.seq, tmp, deparse.level = 0)

    # complexity measure t (t=0 means the graphs is simple)
    mz[g] <- prod(factorial(tmp))
    tz[g] <- m1[g] + log2(mz[g])
    t[g] <- 2 ^ (tz[g])

    # possible pairwise shifts in each edge sequence
    if (m1[g] == m) {
      shifts[g] = 1
    } else {
      shifts[g] = 2 ^ (m - m1[g])
    }

    # total number of graphs given edge sequence
    tot[g] <-  shifts[g] * edge.perm[g]

    # find the simple graphs
    if (m1[g] == 0 & m3[g] == 0) {
      simple[g] = 1
    } else {
      simple[g] = 0
    }
  }

  #total number of graphs given all edge sequences
  Ntot <-
    sum(tot)  # should be equal to: factorial(2*m)/prod(factorial(degvec))

  # probability of each multigraph under RSM
  prob.rsm <- vector()
  for (i in 1:nrow(edge.seq)) {
    tmp <- factorial(m) / prod(factorial(m.seq[i,]))
    prob.rsm[i] <- ((2 ^ (m2[i])) * tmp) / Ntot
  }
  # expected values of some statistics using probability distribution prob.rsm
  lg.prob <- -log2(prob.rsm)
  ElgP <- sum(prob.rsm * lg.prob) # currently not in output
  Em1 <- sum(prob.rsm * m1)
  Em2 <- sum(prob.rsm * m2)
  Vm1 <- sum(prob.rsm * m1 ^ 2) - Em1 ^ 2
  Vm2 <- sum(prob.rsm * m2 ^ 2) - Em2 ^ 2
  Em.seq <- sum(prob.rsm * mz) # currently not in output
  Et <- sum(prob.rsm * tz) # currently not in output

  m.seq <- as.data.frame(m.seq)
  idx <- vector()
  for (i in 1:n) {
    for (j in 1:n) {
      if (i <= j) {
        idx <- c(idx, paste(i, j, sep = ""))
      }
    }
  }
  colnames(m.seq) <- sprintf("M%s", idx[1:r])

  prob.dists <-
    as.data.frame(cbind(
      'prob.rsm' = prob.rsm,
      'loops' = m1,
      'multiedges' = m2,
      'simple' = simple
    ))

  EM <- cbind(Em1, Em2)
  VarM <- cbind(Vm1, Vm2)
  lower <- EM - 2 * sqrt(VarM)
  upper <- EM + 2 * sqrt(VarM)

  out.M <- as.data.frame(rbind(EM, VarM, upper, lower))
  colnames(out.M) <- NULL
  rownames(out.M) <-
    c("Expected", "Variance", "Upper 95%", "Lower 95%")
  colnames(out.M) <-  c("M1", "M2")
  out.M <- round(out.M, 3)

  rsm.out <-
    list('m.seq' = m.seq,
         'prob.dists' = prob.dists,
         'M' = out.M)

  return(rsm.out)
}
