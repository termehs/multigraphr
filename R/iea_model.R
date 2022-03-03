#' @title Independent edge assignment model for multigraphs
#' @description Summary of estimated statistics for analyzing
#' global structure of random multigraphs under the independent edge assignment model
#' given observed adjacency matrix.
#'
#' Two versions of the IEA model are implemented,
#' both of which can be used to approximate the RSM model:
#' \cr
#'   1. independent edge assignment of stubs (IEAS) given an edge probability sequence\cr
#'   2. independent stub assignment (ISA) given a stub probability sequence \cr
#' @param adj matrix of integers representing graph adjacency matrix.
#' @param type equals \code{'graph'} if adjacency matrix is for graphs (default),
#' equals \code{'multigraph'} if it is the equivalence of the adjacency matrix for multigraphs
#' (with matrix diagonal representing loops double counted).
#' @param model character string representing which IEA model: either \code{'IEAS'} (default) or \code{'ISA'}.
#' @param K  upper limit for the complexity statistics \emph{R(k)},
#' \emph{k}=(0,1,...,\code{K}), representing the sequence of frequencies of edge multiplicities
#' (default is maximum observed in adjacency matrix).
#' @param apx logical (default = \code{'FALSE'}). if \code{'TRUE'}, the IEA model is used to approximate
#' the statistics under the random stub matching model given observed degree sequence.
#' @param p.seq if \code{model = ISA} and \code{apx = FALSE}, then specify this numerical vector of
#' stub assignment probabilities.
#' @return
#' \item{nr.multigraphs}{Number of unique multigraphs possible.}
#' \item{M}{Summary and interval estimates for \emph{number of loops} (\code{M1}) and
#' \emph{number of multiple edges} (\code{M2}).}
#' \item{R}{Summary and interval estimates for frequencies of edge multiplicities
#' \code{R1},\code{R2},...,\code{RK}, where \code{K} is a function argument.}
#' @details When using the IEAS model: \cr If the IEAS model is used
#' as an approximation to the RSM model, then the edge assignment probabilities are estimated
#' by using the observed degree sequence. Otherwise, the edge assignment probabilities are
#' estimated by using the observed edge multiplicities  (maximum likelihood estimates). \cr
#'
#' When using the ISA model: \cr If the ISA model is used
#' as an approximation to the RSM model, then the stub assignment probabilities are estimated by using
#' the observed degree sequence over \emph{2m}.  Otherwise, a sequence containing the stub assignment
#' probabilities (for example based on prior belief) should be given as argument \code{p.seq}.
#' @author Termeh Shafie
#' @seealso \code{\link{get_degree_seq}}, \code{\link{get_edge_multip_seq}}, \code{\link{iea_model}}
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#'
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs.
#' \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#'\cr
#'
#' Shafie, T., Schoch, D. (2021). Multiplexity analysis of networks using multigraph representations.
#' \emph{Statistical Methods & Applications} 30, 1425–1444.


#' @examples
#' # Adjacency matrix of a small graph on 3 nodes
#' A <-  matrix(c(1, 1, 0,
#'                1, 2, 2,
#'                0, 2, 0),
#'              nrow = 3, ncol = 3)
#'
#' # When the IEAS model is used
#' iea_model(adj = A , type = 'graph', model = 'IEAS', K = 0, apx = FALSE)
#'
#' # When the IEAS model is used to approximate the RSM model
#' iea_model(adj = A , type = 'graph', model = 'IEAS', K = 0, apx = TRUE)
#'
#' # When the ISA model is used to approximate the RSM model
#' iea_model(adj = A , type = 'graph', model = 'ISA', K = 0, apx = TRUE)
#'
#' # When the ISA model is used with a pre-specified stub assignment probabilities
#'iea_model(adj = A , type = 'graph', model = 'ISA', K = 0, apx = FALSE, p.seq = c(1/3, 1/3, 1/3))
#' @export
#'
iea_model <- function(adj, type = 'multigraph',  model = 'IEAS', K = 0, apx = FALSE, p.seq = NULL) {
  n <- dim(adj)[1]
  r <- choose(n + 1, 2)

  if (type == 'multigraph') {
    if (sum(diag(adj)) %% 2 == 1)
      stop("not an adjacency matrix for multigraphs
           with diagonal elements double counted,
           consider type 'graph' instead.")
    if (sum(adj) %% 2 == 1)
      stop("sum of adjacency matrix must be even")
    m.mat <- adj - 0.5 * diag(diag(adj))
    m.mat[lower.tri(m.mat, diag = FALSE)] <- 0
    m.seq <-  m.mat[upper.tri(m.mat, diag = TRUE)]
    m <- sum(m.seq)
  } else if (type == 'graph') {
    m.mat <- adj
    m.mat[lower.tri(m.mat, diag = FALSE)] <- 0
    m.seq <-  m.mat[upper.tri(m.mat, diag = TRUE)]
    m <- sum(m.seq)
  }

  if (K == 0) {
    K <- max(m.mat)
  } else if (K > m) {
    stop(paste0("K exceeds number of edges = ", m))
  }
  else{
    K <- K
  }

  # possible number of outcomes for edge multiplicity sequences
  mg.outcomes <- choose(m + r - 1, r - 1)

  if (model == 'IEAS'){
    if (!apx) {
      # edge assignment probabilities (Q) as a  matrix and a sequence
      Q.mat <- m.mat / m
      Q.seq <- m.seq / m
    } else if (apx) {
      deg.seq <- get_degree_seq(adj, type)
      Q.mat <- matrix(0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j) {
            Q.mat[i, j] <- deg.seq[i] * (deg.seq[i] - 1) / (2 * m * (2 * m - 1))
          } else if (i != j) {
            Q.mat[i, j] <- 2 * deg.seq[i] * (deg.seq[j]) / (2 * m * (2 * m - 1))
          } else {
            Q.mat[i, j] <- 0
          }
        }
      }
      Q.seq <-  t(Q.mat)[lower.tri(Q.mat, diag = TRUE)]
    }
  } else if (model == 'ISA'){
    if (!apx) {
      if(is.null(p.seq)){
        stop("p.seq must be specified when apx = FALSE")
      }
      Q.mat <- matrix(0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j) {
            Q.mat[i, j] <- p.seq[i]^2
          } else if (i != j) {
            Q.mat[i, j] <- 2*p.seq[i]*p.seq[j]
          } else {
            Q.mat[i, j] <- 0
          }
        }
      }
      Q.seq <-  t(Q.mat)[lower.tri(Q.mat, diag = TRUE)]
    } else if (apx) {
      deg.seq <- get_degree_seq(adj, type)
      p.seq <- deg.seq/(2*m)
      Q.mat <- matrix(0, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j) {
            Q.mat[i, j] <- p.seq[i]^2
          } else if (i != j) {
            Q.mat[i, j] <- 2*p.seq[i]*p.seq[j]
          } else {
            Q.mat[i, j] <- 0
          }
        }
      }
      Q.seq <-  t(Q.mat)[lower.tri(Q.mat, diag = TRUE)]
    }
  } else{
    stop("either the IEAS or ISA model must be specfied")
  }

  # complexity and simplicity statistics under IEA (using derived formulas)
  # m1 = number of loops
  # m2 = number of non-loops
  m1 <- sum(diag(m.mat))
  m2 <- sum(m.mat) - sum(diag(m.mat))

  Em1 <- m * sum(diag(Q.mat))
  Q.upmat <- Q.mat
  Q.upmat[lower.tri(Q.upmat, diag = TRUE)] <- 0
  Em2 <- m * sum(Q.upmat)
  # alt Em2=m-E(m1)

  cov2 <- vector()
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        cov1 <- Q.mat[i, i] * Q.mat[j, j]
        cov2 <- c(cov2, cov1)
      }
    }
  }

  Vm1 <-  m * (sum(diag(Q.mat) * (1 - diag(Q.mat))) - sum(cov2))
  Vm2 <- Vm1

  obs <- cbind(m1, m2)
  EM <- cbind(Em1, Em2)
  VarM <- cbind(Vm1, Vm2)
  lower <- EM - 2 * sqrt(VarM)
  upper <- EM + 2 * sqrt(VarM)

  out.M <- as.data.frame(rbind(obs, EM, VarM, upper, lower))
  colnames(out.M) <- NULL
  rownames(out.M) <-
    c("Observed", "Expected", "Variance", "Upper 95%", "Lower 95%")
  colnames(out.M) <-  c("M1", "M2")
  out.M <- round(out.M, 3)

  # Rk = frequencies of sites with multiplicities k
  R <- vector()
  for (k in 0:K) {
    R[k + 1] <- sum(m.seq == k)
  }

  ER <- vector()
  for (k in 0:K) {
    ER[k + 1] <- sum(choose(m, k) * Q.seq ^ k * (1 - Q.seq) ^ (m - k))
  }

  VarR <- vector()
  lower95 <- vector()
  upper95 <- vector()
  for (k in 0:K) {
    covRk <- matrix(0, r, r)
    for (i in 1:r) {
      for (j in 1:r)
        if (i != j) {
          covRk[i, j] <-
            (choose(m, k) * choose(m - k, k)) * (Q.seq[i] ^ k) * (Q.seq[j] ^ k) * (1 -
                                                                                     Q.seq[i] - Q.seq[j]) ^ (m - 2 * k) -
            ((
              choose(m, k) * (Q.seq[i] ^ k) * (1 - Q.seq[i]) ^ (m - k) * (choose(m, k) *
                                                                            (Q.seq[j] ^ k) * (1 - Q.seq[j]) ^ (m - k))
            ))
        } else{
          covRk[i, j] <- (choose(m, k) * (Q.seq[i] ^ k) * (1 - Q.seq[i]) ^ (m - k)) *
            (1 - (choose(m, k) * Q.seq[i] ^ k * ((1 - Q.seq[i]) ^ (m - k))))
        }
    }
    VarR[k + 1] <- sum(covRk)
    lower95[k + 1] <- ER[k + 1] - 2 * sqrt(VarR[k + 1])
    upper95[k + 1] <- ER[k + 1] + 2 * sqrt(VarR[k + 1])
  }

  out.R <-
    as.data.frame(round(rbind(R, ER, VarR, upper95, lower95), 3))
  colnames(out.R) <- sprintf("R%d", 0:K)
  rownames(out.R) <-
    c("Observed", "Expected", "Variance", "Upper 95%", "Lower 95%")


  # the covariance matrix for local edge multiplicities
  sigma <- matrix(0, r, r)
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
        sigma[i, j] <- Q.seq[i] * (1 - Q.seq[j])
      } else {
        sigma[i, j] <- -(Q.seq[i] * Q.seq[j])
      }
    }
  }

  output <- list("nr.multigraphs" = mg.outcomes, "M" = out.M, "R" = out.R)
  return(output)
}
