# expectedResiduals - Excpected Squared Residuals and Beta Moments

Computes expected residuals and beta2 moments under the variational
posterior

## Usage

``` r
expectedResidualSq(B, i, y, mu, Sigma, p, iter, psi)
```

## Arguments

- B:

  List of basis matrices.

- i:

  Curve index.

- y:

  List of observations.

- mu:

  Matrix of posterior beta means.

- Sigma:

  Array of posterior beta covariance matrices.

- p:

  Matrix of inclusion probabilities.

- iter:

  Current iteration index.

- psi:

  Correlation matrix.

## Value

Numeric expected squared values.
