
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Installation

<!-- You can install the released version of multigraphr from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("multigraphr") -->

<!-- ``` -->

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("termehs/multigraphr")
```

# Overview: multigraphs and applicability

Multigraphs are network representations in which multiple edges and edge
loops (self edges) are permitted. These data structures can be either
directly observed or aggregated by classifying or cross-classifying node
attributes into meta nodes. For the latter case, within group edges
correspond to self-edges (for more details see Shafie 2015;2016). See
example below where the original graph with 15 nodes and 12 edges (left)
is aggregated into a small multigraph with 4 nodes corresponding the
available node attributes (right).

![](mg_ex.png)

# Package overview: `multigraphr`

Two probability models for generating undirected random multigraphs are
implemented together with several statistics under these two models.
Multigraphs are represented by their edge multiplicity sequence, where
the edge multiplicity denotes the number of edges at possible vertex
pair sites ordered according to *(1,1) \< (1,2) \<···\< (1,n) \< (2,2)
\< (2,3) \<···\< (n,n)*, where *n* is number of nodes.

Note that some of the functions are only practical for small scale
multigraphs.

## Random multigraph models

The first model is obtained by random stub matching (RSM) given observed
degree sequence of a multigraphs, so that edge assignments to vertex
pair sites are dependent. The second is obtained by independent edge
assignments (IEA) according to a common probability distribution. There
are two ways in which an approximate IEA model can be obtained from an
RSM model, thus facilitating the analysis. These two ways are
independent stub assignment (ISA) and independent edge assignment of
stubs (IEAS) (Shafie 2015;2016).

### Example

``` r
library('multigraphr')
```

Consider a small multigraph example with 3 nodes and the following
adjacency matrix:

``` r
A <-  matrix(c(1, 1, 0, 
               1, 2, 2, 
               0, 2, 0), 
             nrow = 3, ncol = 3)
A
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
#> [2,]    1    2    2
#> [3,]    0    2    0
```

The degree sequence of the multigraph has double counted diagonals
(stubs for loops) and is given by

``` r
D <- get_degree_seq(adj = A, type = 'graph')
D
#> [1] 2 3 7
```

so that number of edges in the multigraph is 6.

The RSM model given observed degree sequence shows there are 7 possible
multigraphs given this fixed degree sequence, as represented by their
multiplicity sequence `m.seq` (each row correspond to the multiplicity
sequence of a unique multigraph):

``` r
rsm_1 <- rsm_model(deg.seq = D)
rsm_1$m.seq
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    1    1    3
#> [2,]    1    0    0    0    3    2
#> [3,]    0    2    0    0    1    3
#> [4,]    0    1    1    1    0    3
#> [5,]    0    1    1    0    2    2
#> [6,]    0    0    2    1    1    2
#> [7,]    0    0    2    0    3    1
```

with probabilities, number of loops, number of multiple edges and
whether they are simple graphs or not:

``` r
rsm_1$prob.dists
#>     prob.rsm loops multiedges simple
#> 1 0.03030303     5          1      0
#> 2 0.06060606     3          3      0
#> 3 0.06060606     3          3      0
#> 4 0.06060606     4          2      0
#> 5 0.36363636     2          4      0
#> 6 0.18181818     3          3      0
#> 7 0.24242424     1          5      0
```

Consider using the IEA model to aproxiamte the RSM model so that edge
assignment probabilities are functions of observed degree sequence. Note
that the outcome space for multigraphs is much bigger than for the RSM
model so the multiplicity sequences are not printed (they can be found
using the function `get_edgemultip_seq` for very small multigraphs and
probabilities can be found using the multinomial distribution). The
following shows the number of multigraphs under the IEA model:

``` r
iea_1 <-   iea_model(adj = A , type = 'graph', K = 0, apx = TRUE)
iea_1$nr.multigraphs
#> [1] 462
```

## Complexity statistics

The statistic are complexity statistics such as number of loops
(indicator of e.g. homophily) and number of multiple edges (indicator of
e.g. multiplexity/interlocking), together with their probability
distributions, moments and interval estimates.

### Example (cont’d)

First two moments of statistics ‘number of loops’ and ‘number of
multiple edges’ under RSM are given by

``` r
rsm_1$stat.moms
#>   E(loops)  V(loops) E(multiedges) V(multiedges)
#> 1 2.272727 0.9862259      3.727273     0.9862259
```

which are calculated using the probability distributions of the
statistics (no closed formulas exist for these moments). Under the IEA
model, moments of these statistics (*M1* = number of loops, *M2* =
number of multiple edges), together with the complexity statistic *R\_k*
representing the sequence of frequencies of edge sites with
multiplicities *0,1,…,k*, are found using derived formulas:

``` r
iea_1$M
#>               M1    M2
#> Observed   3.000 3.000
#> Expected   2.273 7.455
#> Variance   1.412 1.412
#> Upper 95%  4.649 9.831
#> Lower 95% -0.104 5.078
iea_1$R
#>              R0     R1     R2
#> Observed  2.000  2.000  2.000
#> Expected  2.674  1.588  1.030
#> Variance  0.575  1.129  0.760
#> Upper 95% 4.191  3.713  2.773
#> Lower 95% 1.156 -0.537 -0.713
```

The inteval estimates can then be visualised as box plots to detect
discrepancies between observed and expected, and to detect overlap and
potential depndencies between different types of edges.

## Goodness of fit tests

Goodness of fits tests of multigraph models using Pearson (S) and
information divergence (A) test statistics under the random stub
matching (RSM) and by independent edge assignments (IEA) model, where
the latter is either independent edge assignments of stubs (IEAS) or
independent stub assignment (ISA). The tests are performed using
goodness-of-fit measures between the edge multiplicity sequence of an
observed multigraph, and the expected multiplicity sequence according to
a simple or composite hypothesis.

Probability distributions of test statistics, summary of tests, moments
of tests statistics. adjusted test statistics, critical values,
significance level according to asymptotic distribution, and power of
tests are given.

## Theoretical Background

For more details regarding the theoretical background of the package,
consult the following literature which the package is based on:

  - Shafie, T. (2015). A multigraph approach to social network analysis.
    *Journal of Social Structure*, 16.
    [Link](https://www.exeley.com/journal_of_social_structure/doi/10.21307/joss-2019-011)
  - Shafie, T. (2016). Analyzing local and global properties of
    multigraphs. *The Journal of Mathematical Sociology*, 40(4),
    239-264.
    [Link](https://www.tandfonline.com/doi/abs/10.1080/0022250X.2016.1219732?journalCode=gmas20)

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- library(multigraphr) -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
