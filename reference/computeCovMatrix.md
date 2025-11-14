# computeCovMatrix - Gaussian Process Covariance Matrix

Generates a covariance matrix using the specified Gaussian process,
scaled by sigma2.

## Usage

``` r
computeCovMatrix(Xt, Xp, sigma = 1, w = 1)
```

## Arguments

- Xt:

  Numeric vector of time points.

- Xp:

  Numeric vector of time points.

- sigma:

  Numeric variance parameter.

- w:

  Numeric decay parameter.

## Value

A covariance matrix of dimension length(Xt) x length(Xp)
