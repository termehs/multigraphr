#' @title Edge multiplicity sequences
#' @description  Finds all unique multigraphs represented  by their edge multiplicity sequences
#' @param adj Integer matrix
#' @param type Equals 'graph' if adjacency matrix is for graphs (default),
#' equals 'multigraph' if it is the equivalence of the adjacency matrix for multigraphs
#' (with the elements of the matrix diagonal double counted)
#' @return All unique edge multiplicity sequences as a dataframe
#' @details  To be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#'

get_edgemultip_seq <- function(adj, type = 'multigraph') {
  n <- dim(adj)[1]
  r <- choose(n + 1, 2)

  if (type == 'multigraph') {
    deg.seq <- sort(rowSums(adj))
  } else if (type == 'graph') {
    deg.seq <- sort(rowSums(adj) + diag(adj))
  }

  m <- sum(deg.seq) / 2
  k <- m - 1

  #initial edge sequence (read as labelled 2-tuples of connected nodes connected)
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
      Wz <- z[i, ]
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
  # ordered edge multiplicity sequence = multips.star

  m.seq <- vector()

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
        b <- rbind(b, as.integer(z[j,] == z[i,])) #for multiple edges
        c <- rbind(c, as.integer(z[j,] == f[i,])) #for loops
      }
      if (length(rowSums(b) == 2) > 1) {
        A[z[j, 1], z[j, 2]] <- length(which(rowSums(b) == 2))
      } else {
        A[z[j, 1], z[j, 2]] <- 1
      }
    }

    # the edge multiplicity sequence multips
    tmp <- t(A)[lower.tri(t(A), diag = TRUE)]
    m.seq <- rbind(m.seq, tmp, deparse.level = 0)
  }
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
  return(m.seq)
}
