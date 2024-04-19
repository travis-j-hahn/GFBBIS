
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GFBBIS

<!-- badges: start -->
<!-- badges: end -->

GFBBIS is an implementation of gradient and gradient free (coming) black
box importance sampling originally from Qiang Liu, & Jason D. Lee.
(2016). Black-box Importance Sampling.
<https://arxiv.org/abs/1610.05247>

## Installation

You can install the development version of GFBBIS like so:

``` r
devtools::install_github("travis-j-hahn/GFBBIS")
```

## Example

Here is the basic syntax for the gradient version of black box
importance sampling:

``` r
library(GFBBIS)
## basic example code
out = BBIS(theta,theta_grads,1000,kernel='rbf')

print(out$adj_mean)
```

A full example can be found in examples/logit_example.R
