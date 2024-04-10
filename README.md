
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GFLglm v0.1.0

<!-- badges: start -->
<!-- badges: end -->

This package gives generalized fused Lasso estimator for generalized
linear models via coordinate descent algorithm. Furthermore, this
package includes the R package
[GFLlogit](https://github.com/ohishim/GFLlogit).

**cite this package of the latest version**:  
Ohishi, M. (2024). GFLglm: Generalized fused Lasso for grouped data in
generalized linear models. R package version 0.1.0.
<https://github.com/ohishim/GFLglm>

## Installation

You can install the current version of `GFLglm` like so:

``` r
devtools::install_github("ohishim/GFLglm")
```

## Example

This example requires the following packages:

``` r
library(GFLglm)
library(magrittr)
```

This package has the dataset `crimetko` used in Ohishi (2024) like this

``` r
head(crimetko)
#>   year area crime pop group
#> 1 2017    1   714   8     1
#> 2 2017    1    76   5     1
#> 3 2017    1    44   4     1
#> 4 2017    1   275  22     1
#> 5 2017    1    61   4     1
#> 6 2017    1    47 526     1
```

The `crime` is the number of crimes of which the original data was
collected by the Metropolitan Police Department, available at TOKYO OPEN
DATA (<https://portal.data.metro.tokyo.lg.jp/>). This variable is
arranged and used the following production:

- Tokyo Metropolitan Government & Metropolitan Police Department. The
  number of recognized cases by region, crime type, and method (yearly
  total; in Japanese),
  <https://creativecommons.org/licenses/by/4.0/deed.en>.

The `pop` is population of which the original data was obtained from the
results of the population census, as provided in e-Stat
(<https://www.e-stat.go.jp/en>). This variable is arranged and used.

This package also has the dataset `adj` of adjacent relationship of
groups in `crimetko` like this

``` r
head(adj)
#>   g adj
#> 1 1   2
#> 2 1   3
#> 3 1   4
#> 4 1   5
#> 5 1   6
#> 6 1  54
```

For example, group 1 is adjacent to groups 2, 3, 4, 5, 6, 54, … For such
a data, Poisson regression with generalized fused Lasso can be executed
as

``` r
y <- crimetko$crime                           # the response variable  
group <- crimetko$group                       # group id
D <- adj %$% split(adj, g)                    # adjacent relationship
q <- crimetko$pop %>% inset(.==0, 1) %>% log  # offset for log-link

res <- GFLglm(y, group, D, "Poisson", offset=q)
```

The fourth parameter specifies the error distribution which can take
“Gaussian”, “Binomial”, “Poisson”, “Negative.Binomial” (or “NB”),
“Gamma”, or “Inverse.Gaussian” (or “IG”). For “Gamma”, “NB”, and “IG”,
there is a option using log-link by setting `link = "log"`.

## Reference

1.  Ohishi, M. (2024). Generalized fused Lasso for grouped data in
    generalized linear models. *Hiroshima Statistical Research Group
    Technical Report*, TR-No. 24-02, Hiroshima University.
    \[[PDF](http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR24/TR24-02.pdf)\]
