
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview: multigraphs and applicability

Multigraphs are network representations in which multiple edges and edge
loops (self edges) are permitted. These data structures can be either
directly observed or aggregated by classifying or cross-classifying node
attributes into meta nodes. For the latter case, within group edges
correspond to self-edges (for more details see Shafie 2015;2016). See
example below where the original graph with four node attributes (left)
are aggregated into a multigraph (right).

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
degree sequence of a multigrpahs, so that edge assignments to vertex
pair sites are dependent. The second is obtained by independent edge
assignments (IEA) according to a common probability distribution. There
are two ways in which an approximate IEA model can be obtained from an
RSM model faciliating the ananlysis. These two ways are indepndent stub
assignment (ISA) and independent edge assignment of stubs (IEAS) (Shafie
2015;2016).

## Complexity statistics

The statistic are complexity statistics such as number of loops
(indicator of homophily) and number of multiple edges (indicator of
multiplexity/interlocking of ties), together with their probability
distributions, moments and interval estimates.

## Goodness of fit tests

Goodness of fits tests of multigraph models using Pearson (S) and
information divergence (A) test statistics under the random stub
matching (RSM) and by independent edge assignments (IEA) model, where
the latter is either independent edge assignments of stubs (IEAS) or
independent stub assignment (ISA). The tests are performed using
goodness-of-fit measures between the edge multiplicity sequence of an
observed multigraph, and the expected multiplicity sequence according to
a simple or composite hypothesis.

Porbability distributions of test statistics, summary of tests, moments
of tests statistics. adjusted test statistics, critical values,
significance level according to asymptotic distribution, and power of
tests are
given.

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

## Example

## Theoretical Background

`multigraphr` is based on the following published papers. For more
details regarding the theoretical background of the package, consult the
following literature:

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
