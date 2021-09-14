#' @title Edge multiplicity sequences of multigraphs given fixed degrees
#' @description  Given a degree sequence, this function finds all unique multigraphs
#' represented by their edge multiplicity sequences.
#' @param deg.seq vector of integers with the sum equal to 2\code{m} representing
#' the degree sequence of the multigraph.
#' @return All unique edge multiplicity sequences as rows in a data frame.
#' Each row in the data frame represents a unique multigraph given the degree sequence.
#' @details   Multigraphs are represented by their edge multiplicity sequence \strong{M} with elements \emph{M(i,j)},
#' denoting edge multiplicity at possible vertex pair sites \emph{(i,j)} ordered according to\cr
#' \emph{(1,1) < (1,2) <···< (1,n) < (2,2) < (2,3) <···< (n,n)}, \cr
#' where \emph{n} is number of nodes.
#'
#' Only practical for small multigraphs.
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @seealso \code{\link{get_degree_seq}}
#' @examples
#' # Adjacency matrix for undirected network with 3 nodes
#'  A <-  matrix(c(0, 1, 2,
#'                 1, 2, 1,
#'                 2, 1, 2), nrow=3, ncol=3)
#' deg <- get_degree_seq(A, 'multigraph')
#' get_edge_multip_seq(deg)
#' @export

get_edge_multip_seq <- function(deg.seq) {
  if (sum(deg.seq) %% 2 == 1)
    stop("sum of degree sequence must be an even number")
  n <- length(deg.seq)
  r <- choose(n + 1, 2)
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

  # initial edge list
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

  # edge multiplicity sequence = m.seq
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

    # the edge multiplicity sequence
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
