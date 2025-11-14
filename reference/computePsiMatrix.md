# computePsiMatrix - Correlation (Covariance) Matrix

Computes Ψ, the correlation matrix, using an exponential kernel.

## Usage

``` r
computePsiMatrix(Xt, Xp, w = 1)
```

## Arguments

- Xt:

  Numeric vector of time points.

- Xp:

  Numeric vector of time points (often identical to Xt).

- w:

  Numeric scalar decay parameter, controlling smoothness.

## Value

A covariance matrix of dimension length(Xt) x length(Xp)
