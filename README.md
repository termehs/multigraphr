
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/networkdata)](https://cran.r-project.org/package=multigrapr)
<!-- badges: end -->

# Package overview: `multigraphr`

## Installation

<!-- You can install the released version of multigraphr from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("multigraphr") -->

<!-- ``` -->

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("termehs/multigraphr")
```

# Multigraphs and applicability

Multigraphs are network representations in which multiple edges and edge
loops (self edges) are permitted. These data structures can be either
directly observed or aggregated by classifying or cross-classifying node
attributes into meta nodes. For the latter case, within group edges
correspond to self-edges (for more details see Shafie 2015;2016). See
example below where the original graph with 15 nodes and 12 edges (left)
is aggregated into a small multigraph with 4 nodes corresponding the
available node attributes (right).

![](mg_ex.png) Multigraphs are represented by their edge multiplicity
sequence, where the edge multiplicity denotes the number of edges at
possible vertex pair sites ordered according to *(1,1) \< (1,2) \<···\<
(1,n) \< (2,2) \< (2,3) \<···\< (n,n)*, where *n* is number of nodes.

Two probability models for generating undirected random multigraphs are
implemented in the package together with several statistics under these
two models. Moreover, functions for goodness of fit tests are available
for the presented models.

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

Consider a small graph on 3 nodes and the following adjacency matrix:

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

Consider using the IEA model to approximate the RSM model so that edge
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

Under the RSM model, the first two moments and interval estimates of the
statistics *M1* = ‘number of loops’ and *M2* = ‘number of multiple
edges’ are given by

``` r
rsm_1$M
#>              M1    M2
#> Expected  2.273 3.727
#> Variance  0.986 0.986
#> Upper 95% 4.259 5.713
#> Lower 95% 0.287 1.741
```

which are calculated using the numerically found probability
distributions under RSM (no analytical solutions exist for these
moments).

Under the IEA model, moments of these statistics, together with the
complexity statistic *R\_k* representing the sequence of frequencies of
edge sites with multiplicities *0,1,…,k*, are found using derived
formulas:

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

The interval estimates can then be visualised to detect discrepancies
between observed and expected, and to detect overlap and potential
dependencies between different types of edges.

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

### Example

Goodness of fit tests for multigraphs with *n=4* nodes and *m=10* edges.

  - Testing a simple IEAS hypothesis with degree sequence *(5,5,5,5)*
    against a RSM model with degree sequence
*(14,2,2,2)*:

<!-- end list -->

``` r
gof_1 <- gof_multigraph(m = 10, model = 'RSM', deg.mod = c(14,2,2,2), hyp = 'IEAS', deg.hyp = c(5,5,5,5))
gof_1
#> $probS
#>      S=s  P(S=s)  P(S<s)
#> 1  29.52 0.34675 0.34675
#> 2  42.82 0.41610 0.76285
#> 3  45.48 0.10402 0.86687
#> 4  61.44 0.06935 0.93622
#> 5  62.58 0.03467 0.97090
#> 6  64.48 0.01734 0.98824
#> 7  65.24 0.00867 0.99690
#> 8  85.38 0.00165 0.99856
#> 9  88.04 0.00124 0.99979
#> 10 88.80 0.00021 1.00000
#> 
#> $probA
#>         A=a  P(A=a)  P(A<a)
#> 1  21.24967 0.34675 0.34675
#> 2  22.54115 0.41610 0.76285
#> 3  27.00793 0.06935 0.83220
#> 4  27.14628 0.10402 0.93622
#> 5  28.84047 0.03467 0.97090
#> 6  32.55310 0.01734 0.98824
#> 7  33.44560 0.00867 0.99690
#> 8  34.58208 0.00165 0.99856
#> 9  39.18721 0.00124 0.99979
#> 10 40.07985 0.00021 1.00000
#> 
#> $summmary
#>   stat      exp       var       cv   alpha stat>cv cv(stat) stat>cv(stat)
#> 1    S 41.16706 115.37420 17.48528 0.04164       1 62.64954        0.0291
#> 2    A 23.41280   8.19897 17.48528 0.04164       1 29.13957        0.0291
#> 
#> $adjusted.stats
#>                   S    Sprim      Sbis         A    Aprim     Abis
#> Exp        41.16706 41.16706  41.16706 23.412802 23.41280  23.4128
#> Var       115.37420 82.66960 376.60594  8.198965 47.66603 121.8132
#> Best Adj.   0.00000  1.00000   0.00000  0.000000  1.00000   0.0000
#> 
#> $adjusted.chi2
#>      df(S) df(A)
#> [1,]    41    23
#> 
#> $power.apx
#>     Sprim Sbis   Aprim Abis
#> 1 0.00046    1 0.19952    1
#> 
#> $degrees.of.freedom
#> [1] 9
```

  - Testing a composite IEAS hypothesis against a RSM model with degree
    sequence
*(14,2,2,2)*:

<!-- end list -->

``` r
gof_2 <- gof_multigraph(m = 10, model = 'RSM', deg.mod = c(14,2,2,2), hyp = 'IEAS', deg.hyp = 0)
gof_2
#> $probS
#>        S=s  P(S=s)  P(S<s)
#> 1  1.48352 0.34675 0.34675
#> 2  4.04121 0.41610 0.76285
#> 3  8.37363 0.06935 0.83220
#> 4 14.48077 0.00165 0.83385
#> 5 19.23077 0.01734 0.85119
#> 6 19.64835 0.10402 0.95521
#> 7 22.62363 0.03467 0.98989
#> 8 38.23077 0.00991 0.99979
#> 9 57.23077 0.00021 1.00000
#> 
#> $probA
#>         A=a  P(A=a)  P(A<a)
#> 1   2.22362 0.34675 0.34675
#> 2   3.21690 0.41610 0.76285
#> 3   7.38562 0.06935 0.83220
#> 4   8.76207 0.10402 0.93622
#> 5  10.15821 0.03467 0.97090
#> 6  12.93080 0.01734 0.98824
#> 7  14.66173 0.00165 0.98989
#> 8  15.70339 0.00867 0.99856
#> 9  20.20690 0.00124 0.99979
#> 10 22.97949 0.00021 1.00000
#> 
#> $summmary
#>   stat     exp      var      cv   alpha stat>cv cv(stat) stat>cv(stat)
#> 1    S 6.35294 51.88277 12.9282 0.04419  0.1678 20.75888       0.04479
#> 2    A 4.29978  8.88590 12.9282 0.04419  0.0291 10.26162       0.02910
#> 
#> $adjusted.stats
#>                   S     Sprim      Sbis        A    Aprim     Abis
#> Exp        6.352943  6.352943  6.352943 4.299776 4.299776 4.299776
#> Var       51.882769 13.453295 13.453295 8.885904 9.244037 6.162691
#> Best Adj.  0.000000  0.000000  0.000000 0.000000 1.000000 0.000000
#> 
#> $adjusted.chi2
#>      df(S) df(A)
#> [1,]     6     4
#> 
#> $power.apx
#>     Sprim    Sbis   Aprim    Abis
#> 1 0.94255 0.99988 0.98285 0.99903
#> 
#> $degrees.of.freedom
#> [1] 6
```

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
