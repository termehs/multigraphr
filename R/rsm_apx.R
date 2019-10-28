#' @title IEA approximation of the statistics under the RSM model
#' @description Summary of statistics estimated for
#' global structure of random multigraphs under an IEA approximation of the RSM model.
#' The degree sequence is used to estimate the edge assignment probabilities.
#' @param deg.seq vector of integers representing the degree sequence of the multigraph
#' @return List including number of multigraphs under the IEA approximation,
#' summary and interval estomates for number of loops and number of non-loos (M_1 and M__2),
#' summary and interval estimates for frequencies of edge multiplicites R_k
#' @details  To be completed
#' @author Termeh Shafie
#' @references Shafie, T. (2015). A Multigraph Approach to Social Network Analysis. *Journal of Social Structure*, 16.
#' Shafie, T. (2016). Analyzing Local and Global Properties of Multigraphs. *The Journal of Mathematical Sociology*, 40(4), 239-264.
#'
#'
rsm_apx <- function(deg.seq){
  length(deg.seq)
  r <- choose(n+1,2)

  # edge assignment probabilities given degree sequence
  Q.mat <- matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {
        Q.mat[i,j] <- deg.seq[i]*(deg.seq[i]-1)/(2*m*(2*m-1))
      } else if (i<j){
        Q.mat[i,j] <- 2*deg.seq[i]*(deg.seq[j])/(2*m*(2*m-1))
      } else {
        Q.mat[i,j] <- 0
      }
    }
  }


  # outcome space for the edge multiplcitiy sequences is too large, relies on formulas or the edge multiplcity sequences from RSM
  # define the edge assignment probabilities Q given degree sequence
  # (for the multinomial distribution where m_ij are binomially distributed with prob Q_ij)
  Q.seq <- Q.mat[upper.tri(Q.mat, diag = TRUE)]
  m.seq <- edge_multip_seq(adj, type)
  # find probabilities of RSM outcomes with multinomial distribution (given multips)
  prob.iea <- vector()
  for (i in 1:nrow(m.seq)) {
    prob.iea[i] <- dmultinom(m.seq[i,], m, Q.seq, log = FALSE)
  }
  # complexity and simplicity statistics under IEA (using derived formulas)
  # m1 = number of loops
  # m2 = number of non-loops
  m1 <- sum(diag(m.mat))
  m2 <- sum(m.mat)-sum(diag(m.mat))

  Em1 <- m*sum(diag(Q.mat))
  Em2 <- m*sum(Q.mat[row(Q.mat) != col(Q.mat)])
  # alt Em2=m-E(m1)

  cov2 <- vector()
  for (i in 1:n) {
    for (j in 1:n) {
      if (i!=j) {
        cov1 <- Q.mat[i,i]*Q.mat[j,j]
        cov2 <- c(cov2, cov1)
      }
    }
  }
  Vm1 <-  m*(sum(diag(Q.mat)*(1-diag(Q.mat)))-sum(cov2))
  Vm2 <- Vm1

  # covariance matrix for m1 and m2
  cov <- matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        cov[i,i] <- Q.mat[i,i]*(1-Q.mat[i,i])
      } else{
        cov[i,j] <- -Q.mat[i,i]*Q.mat[j,j]
      }
    }
  }
  Vm1 <- m*sum(cov)
  Vm2 <- Vm1

  obs <- cbind(m1, m2)
  EM <- cbind(Em1, Em2)
  VarM <- cbind(Vm1, Vm2)
  lower <- EM - 2*sqrt(VarM)
  upper <- EM + 2*sqrt(VarM)

  out.M <- as.data.frame(rbind(obs, EM, VarM, upper, lower))
  colnames(out.M) <- NULL
  rownames(out.M) <- c("Observed", "Expected", "Variance", "Upper 95%", "Lower 95%")
  colnames(out.M) <-  c("M1", "M2")
  out.M <- round(out.M,3)

  # Rk = frequencies of sites with multiplicities k
  R <- vector()
  for (k in 0:K) {
    R[k+1] <- sum(m.seq==k)
  }


  ER <- vector()
  for (k in 0:K) {
    ER[k+1] <- sum(choose(m,k)*Q.seq^k*(1-Q.seq)^(m-k))
  }

  VarR <- vector()
  lower95 <- vector()
  upper95 <- vector()
  for (k in 0:K) {
    covRk <- matrix(0,r,r)
    for (i in 1:r) {
      for (j in 1:r)
        if(i != j){
          covRk[i,j] <- (choose(m,k)*choose(m-k,k))*(Q.seq[i]^k)*(Q.seq[j]^k)*(1-Q.seq[i]-Q.seq[j])^(m-2*k)-
            ((choose(m,k)*(Q.seq[i]^k)*(1-Q.seq[i])^(m-k)*(choose(m,k)*(Q.seq[j]^k)*(1-Q.seq[j])^(m-k))))
        } else{
          covRk[i,j] <- (choose(m,k)*(Q.seq[i]^k)*(1-Q.seq[i])^(m-k))*
            (1-(choose(m,k)*Q.seq[i]^k*((1-Q.seq[i])^(m-k))))
        }
    }
    VarR[k+1] <- sum(covRk)
    lower95[k+1] <- ER[k+1] - 2*sqrt(VarR[k+1])
    upper95[k+1] <- ER[k+1] + 2*sqrt(VarR[k+1])
  }

  out.R <- as.data.frame(round(rbind(R,ER,VarR, upper95, lower95),3))
  colnames(out.R) <- sprintf("R%d", 0:K)
  rownames(out.R) <- c("Observed", "Expected", "Variance", "Upper 95%", "Lower 95%")


  # the covariance matrix for local edge multiplicities (delete??)
  sigma <- matrix(0,r,r)
  for (i in 1:r) {
    for (j in 1:r) {
      if (i == j) {
        sigma[i,j] <- Q.seq[i]*(1-Q.seq[j])
      } else {
        sigma[i,j] <- -(Q.seq[i]*Q.seq[j])
      }
    }
  }



}

