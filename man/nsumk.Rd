% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsumk.R
\name{nsumk}
\alias{nsumk}
\title{Ordered n-tuples of non-negative integers summing to k}
\usage{
nsumk(n, k)
}
\arguments{
\item{n}{A positive integer}

\item{k}{A positive integer}
}
\value{
A matrix with nchoosek(k+n-1,n-1) rows and n columns. Each row comprises non-negative integers summing to k.
}
\description{
Finds ordered n-tuples of non-integers summing to k.
Only practical for n < 15.
}
\details{
Useful for finding all possible degree sequences for a network with n nodes
and k/2 number of edges.
}
\examples{
## All possible degree sequences for a network with 4 nodes and 5 edges
D <- nsumk(4, 10)
# Remove isolated nodes
D <- D[-which(rowSums(D == 0) > 0), ]
}
\author{
Termeh Shafie
}
