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
#' \cr
#' #' Shafie, T., Schoch, D. (2021). Multiplexity analysis of networks using multigraph representations.
#' \emph{Statistical Methods & Applications} 30, 1425–1444.
#'
#' @seealso \code{\link{get_degree_seq}}
#' @examples
#' # Adjacency matrix for undirected network with 3 nodes
#' A <- matrix(c(
#'     0, 1, 2,
#'     1, 2, 1,
#'     2, 1, 2
#' ), nrow = 3, ncol = 3)
#' deg <- get_degree_seq(A, "multigraph")
#' get_edge_multip_seq(deg)
#' @export

get_edge_multip_seq <- function(deg.seq) {
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
            # edge.seq <- rbind(edge.seq, s, deparse.level = 0)
            for (i in 1:m) {
                z[i, 1] <- s[2 * i - 1]
                z[i, 2] <- s[2 * i]
            }
            k <- m - 1
        } else {
            k <- k - 1
        }
    }
    edge.seq <- do.call(rbind, edge.seqs)
    # Pre-allocate matrices and vectors
    A <- matrix(0, n, n) # Adjacency matrix for each configuration
    m.seq <- matrix(0, nrow(edge.seq), choose(n, 2) + n) # Pre-allocate for all configurations

    # Iterate over each graph configuration
    for (g in seq_len(nrow(edge.seq))) {
        # Directly use the matrix 'z' constructed from 'edge.seq'
        z <- matrix(edge.seq[g, ], ncol = 2, byrow = TRUE)

        # Reset the adjacency matrix for each configuration
        A[, ] <- 0

        # Vectorized operation to fill the adjacency matrix
        for (i in seq_len(nrow(z))) {
            A[z[i, 1], z[i, 2]] <- A[z[i, 1], z[i, 2]] + 1
        }

        # Extract the lower triangle including the diagonal as the edge multiplicity sequence
        m.seq[g, ] <- t(A)[lower.tri(t(A), diag = TRUE)]
    }
    m.seq <- as.data.frame(m.seq)

    comb <- expand.grid(1:n, 1:n)
    filtered_comb <- comb[comb$Var1 <= comb$Var2, ]
    idx <- paste(filtered_comb$Var1, filtered_comb$Var2, sep = "")

    colnames(m.seq) <- sprintf("M%s", idx[1:r])
    return(m.seq)
}
