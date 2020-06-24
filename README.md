
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview: multigraphs and applicability

Multigraphs are network representations in which multiple edges and edge
loops (self edges) are permitted. These data structures can be either
directly observed or aggregated by classifying or cross-classifying node
attributes into meta nodes. For the latter case, within group edges
correspond to self-edges (for more details see Shafie 2015;2016). See
example below where the original graph with four node attributes (left)
are aggregated into a multigraph
(right).

![](mg_ex.png)

## Package overview: `multigraphr`

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

<!-- ## Example -->

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
