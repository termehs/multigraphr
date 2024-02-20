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
#' D <- c(2, 2, 3, 3) # degree sequence
#' mod1 <- rsm_model(D)
#' mod1$m.seq
#' mod1$prob.dists
#' mod1$M
#'
#' # Given an observed graph
#' A <- matrix(c(
#'     0, 1, 2,
#'     1, 2, 1,
#'     2, 1, 2
#' ), nrow = 3, ncol = 3)
#' D <- get_degree_seq(adj = A, type = "graph")
#' mod2 <- rsm_model(D)
#' @export
#'
rsm_model <- function(deg.seq) {
    if (sum(deg.seq) %% 2 == 1) {
        stop("sum of degree sequence must be an even number")
    }
    n <- length(deg.seq)
    r <- choose(n + 1, 2)
    m <- sum(deg.seq) / 2
    k <- m - 1

    # initial edge sequence (read as labelled 2-tuples of connected nodes connected)
    s <- edge.seq <- rep(seq_len(n), deg.seq)
    edge.seqs <- list(edge.seq)
    # initial edge list
    z <- matrix(edge.seq, ncol = 2, byrow = TRUE)
    # create all possible sequences of edges given degree sequence
    while (k > 0) {
        tz <- t(z)
        W <- c(tz[, k:m])

        tmp <- c(W[1], W[2] + 1)
        W <- sort(W)
        if (length(which(tmp[1] == W)) <= length(which(tmp[2] <= W))) {
            i1 <- which(tmp[1] == W)
            i2 <- which(tmp[2] <= W)
            ni <- seq_len(min(length(i1), length(i2)))
            del <- c(rbind(i1[ni], i2[ni]))
            W1 <- c(rbind(W[i1[ni]], W[i2[ni]]))

            W <- W[-del]
            ss <- c(W1, W)
            s[(2 * m - length(ss) + 1):(2 * m)] <- ss
            edge.seqs[[length(edge.seqs) + 1]] <- s
            for (i in 1:m) {
                z[i, 1] <- s[2 * i - 1]
                z[i, 2] <- s[2 * i]
            }
            k <- m - 1
        } else {
            k <- k - 1
        }
    }
    edge.seq <- do.call(rbind,edge.seqs)
    # for each possible edge sequence/multigraph given degree sequence, count number of loops and number of multiple edges using:
    # edge multiplicity sequence = m.seq
    # ordered edge multiplicity sequence = m.seq.star
    m.seq <- vector()
    m1 <- numeric(nrow(edge.seq))
    m2 <- numeric(nrow(edge.seq))
    m3 <- numeric(nrow(edge.seq))
    edge.perm <- vector()
    mz <- numeric(nrow(edge.seq))
    tz <- numeric(nrow(edge.seq))
    t <- numeric(nrow(edge.seq))
    shifts <- numeric(nrow(edge.seq))
    tot <- numeric(nrow(edge.seq))
    simple <- numeric(nrow(edge.seq))
    m.seq <- matrix(0, nrow(edge.seq), choose(n, 2) + n) # Pre-allocate for all configurations

    for (g in seq_len(nrow(edge.seq))) {
        z <- matrix(edge.seq[g, ], ncol = 2, byrow = TRUE)
        A <- matrix(0, n, n) # adjacency matrix
        for (i in seq_len(nrow(z))) {
            A[z[i, 1], z[i, 2]] <- A[z[i, 1], z[i, 2]] + 1
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
            edge.perm[g] <- factorial(m)
        }

        # the edge multiplicity sequence m.seq
        m.seq[g, ] <- t(A)[lower.tri(t(A), diag = TRUE)]

        # complexity measure t (t=0 means the graphs is simple)
        mz[g] <- prod(factorial(tmp))
        tz[g] <- m1[g] + log2(mz[g])
        t[g] <- 2^(tz[g])

        # possible pairwise shifts in each edge sequence
        if (m1[g] == m) {
            shifts[g] <- 1
        } else {
            shifts[g] <- 2^(m - m1[g])
        }

        # total number of graphs given edge sequence
        tot[g] <- shifts[g] * edge.perm[g]

        # find the simple graphs
        if (m1[g] == 0 && m3[g] == 0) {
            simple[g] <- 1
        } else {
            simple[g] <- 0
        }
    }

    # total number of graphs given all edge sequences
    Ntot <- sum(tot) # should be equal to: factorial(2*m)/prod(factorial(degvec))

    # probability of each multigraph under RSM
    prob.rsm <- vector()
    for (i in seq_len(nrow(edge.seq))) {
        tmp <- factorial(m) / prod(factorial(m.seq[i, ]))
        prob.rsm[i] <- ((2^(m2[i])) * tmp) / Ntot
    }
    # expected values of some statistics using probability distribution prob.rsm
    # lg.prob <- -log2(prob.rsm)
    # ElgP <- sum(prob.rsm * lg.prob) # currently not in output
    Em1 <- sum(prob.rsm * m1)
    Em2 <- sum(prob.rsm * m2)
    Vm1 <- sum(prob.rsm * m1^2) - Em1^2
    Vm2 <- sum(prob.rsm * m2^2) - Em2^2
    # Em.seq <- sum(prob.rsm * mz) # currently not in output
    # Et <- sum(prob.rsm * tz) # currently not in output

    m.seq <- as.data.frame(m.seq)
    comb <- expand.grid(1:n, 1:n)
    filtered_comb <- comb[comb$Var1 <= comb$Var2, ]
    idx <- paste(filtered_comb$Var1, filtered_comb$Var2, sep = "")
    colnames(m.seq) <- sprintf("M%s", idx[1:r])

    prob.dists <-
        as.data.frame(cbind(
            "prob.rsm" = prob.rsm,
            "loops" = m1,
            "multiedges" = m2,
            "simple" = simple
        ))

    EM <- cbind(Em1, Em2)
    VarM <- cbind(Vm1, Vm2)
    lower <- EM - 2 * sqrt(VarM)
    upper <- EM + 2 * sqrt(VarM)

    out.M <- as.data.frame(rbind(EM, VarM, upper, lower))
    colnames(out.M) <- NULL
    rownames(out.M) <-
        c("Expected", "Variance", "Upper 95%", "Lower 95%")
    colnames(out.M) <- c("M1", "M2")
    out.M <- round(out.M, 3)

    rsm.out <-
        list(
            "m.seq" = m.seq,
            "prob.dists" = prob.dists,
            "M" = out.M
        )

    return(rsm.out)
}
