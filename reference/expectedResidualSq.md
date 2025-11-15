# helpers - beta: beta-related expectation helpers

Expected squared residuals and beta^2 moments under the variational
posterior. These are used in the updates of q(sigma^2), q(tau^2), and
q(Z_ki).

## Usage

``` r
expectedResidualSq(B, i, y, mu, Sigma, prob, iter, psi)
```

## Arguments

- B:

  List of basis matrices (B\[i\] is n_i x K).

- i:

  Curve index.

- y:

  List of observations (y\[i\] is length n_i vector).

- mu:

  Matrix of posterior beta means (rows = iterations, cols = beta_k_i).

- Sigma:

  Array of posterior beta covariance matrices (K x K x m).

- prob:

  Matrix of inclusion probabilities (rows = iterations, cols =
  beta_k_i).

- iter:

  Current iteration index.

- psi:

  Correlation matrix for the errors.

## Value

Numeric expected squared (or summed squared) values.
