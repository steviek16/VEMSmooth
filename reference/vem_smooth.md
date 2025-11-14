# vem_smooth - Variational EM Algorithm for Bayesian Basis Selection

Runs the Variational EM (VEM) algorithm for smoothing functional curves
via basis function selection, accounting for correlated errors.

## Usage

``` r
vem_smooth(
  y,
  B,
  Xt = seq(0, 1, length.out = max(length(y))),
  m = length(y),
  K = ncol(B[[1]]),
  mu_ki = 0.5,
  lambda_1 = 1e-10,
  lambda_2 = 1e-10,
  delta_1 = 1e-10,
  delta_2 = 1e-10,
  maxIter = 1000,
  initial_values,
  convergence_threshold = 0.01,
  lower_opt = 0.1
)
```

## Arguments

- y:

  List of length m, each element a numeric vector of observed curve
  values (ie. outcome variable).

- B:

  List of length m, each element of a n_i x K matrix of basis
  evaluations (eg. B-spline values).

- Xt:

  Vector of time points (assumed to be identical for all curves).

- m:

  Number of curves.

- K:

  Number of basis functions.

- mu_ki:

  Hyperparameter for the Beta prior on θ (inclusion probability).

- lambda_1, lambda_2:

  Hyperparameter for the inverse-gamma prior τ².

- delta_1, delta_2:

  Hyperparameters for the inverse-gamma prior on σ².

- maxIter:

  Maximum number of VEM iterations.

- initial_values:

  List of starting values for for w, p, λ₂, and δ₂.

- convergence_threshold:

  ELBO convergence tolerance.

- lower_opt:

  Lower bound for w optimization.

## Value

A list containing fitted variational parameters, ELBO trajectory,
convergence information and final w.

## Details

The algorithm alternates between:

- **E-step:** updating the variartional distributions for β, σ², τ², Z,
  and θ.

- **M-step:** maximizing the ELBO with respect to the correlation
  parameter w.
