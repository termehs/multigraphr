#' @title Complexity Statistics under the IEA model for multigraphs
#' @description Summary of statistics estimated for analysing
#' global structure of random multigraphs under independent edge assignment model
#' given observed adjacency matrix.
#' The edge assignment probabilities are estimated using the observed edge multiplicities
#' (maximum likelihood estimates)
#' @param adj Matrix of integers.
#' @param type Equals 'graph' if adjacency matrix is for graphs (default),
#' equals 'multigraph' if it is the equivalence of the adjacency matrix for multigraphs
#' (with the matrix diagonal double counted).
#' @param K  Upper limit for k in the complexity statistics \emph{R_k} representing the sequence of
#' frequencies of edge sites with multiplicities \emph{0,1,...,k}. Default is maximum observed in adjacency matrix.
#' @param apx logical (default = 'FALSE'). if 'TRUE', the IEA model is used to approximate
#' the statistics under the random stub matching model given observed degree sequence (use function 'get_degree_seq').
#' @return
#' \item{nr.multigraphs}{Number of unique multigraphs possible.}
#' \item{M}{Summary and interval estimates for 'number of loops' and 'number of multiple edges' (\emph{M1} and \emph{M2})).}
#' \item{R}{Summary and interval estimates for frequencies of edge multiplicities \emph{R_k}.}
#' @details  To be completed
#' @author Termeh Shafie
#' @seealso [get_degree_seq]
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. \emph{Journal of Social Structure}, 16.
#' \cr
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. \emph{The Journal of Mathematical Sociology}, 40(4), 239-264.
#' @examples
#' ## Adjacency matrix of a small graph on 3 nodes
#' A <-  matrix(c(1, 1, 0,
#'                1, 2, 2,
#'                0, 2, 0),
#'              nrow = 3, ncol = 3)
#'iea_model(adj = A , type = 'graph', K = 0, apx = FALSE)
#' @export
#'
  iea_model <- function(adj, type = 'multigraph' ,  K = 0, apx = FALSE) {
    n <- dim(adj)[1]
    r <- choose(n + 1, 2)

    if (type == 'multigraph') {
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

    if (apx == FALSE) {
      # edge assignment probabilities (Q) as a  matrix and a sequence
      Q.mat <- m.mat / m
      Q.seq <- m.seq / m
    } else if (apx == TRUE) {
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

    # complexity and simplicity statistics under IEA (using derived formulas)
    # m1 = number of loops
    # m2 = number of non-loops
    m1 <- sum(diag(m.mat))
    m2 <- sum(m.mat) - sum(diag(m.mat))

    Em1 <- m * sum(diag(Q.mat))
    Em2 <- m * sum(Q.mat[row(Q.mat) != col(Q.mat)])
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
