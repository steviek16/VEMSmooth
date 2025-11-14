# expectedLogLikelihood

Internal helper functions for computing the Evidence Lower Bound (ELBO)
components within the Variational EM algorithm for correlated functional
data.

Each function corresponds to a specific ELBO sub-term:

- expectedLogLikelihood(): expected log-likelihood for each curve.

- elboInclusionTerm(): contribution from inclusion indicators (Z terms).

- elboBetaTerm(): contribution from β coefficients.

- elboThetaTerm(): contribution from θ prior parameters.

- elboSigmaTerm(): contribution from σ² prior.

- elboTauTerm(): contribution from τ² prior.

All functions are internal and are not meant to be called directly by
users.

Computes the expected log-likelihood term for curve i given correlated
errors.

## Usage

``` r
expectedLogLikelihood(
  y,
  ni,
  B,
  i,
  iter,
  delta1_q,
  delta2_values,
  mu_beta_values,
  Sigma_beta,
  p_values,
  psi
)
```

## Arguments

- y:

  List of observed data for each curve.

- ni:

  Vector of sample sizes per curve.

- B:

  List of basis matrices.

- i:

  Curve index.

- iter:

  Current iteration index.

- delta1_q:

  Shape parameter of sigma²’s variational posterior.

- delta2_values:

  Rate parameters of sigma²’s variational posterior.

- mu_beta_values:

  Posterior mean matrix of beta.

- Sigma_beta:

  Posterior covariance array of beta.

- p_values:

  Matrix of inclusion probabilities.

- psi:

  Correlation matrix (Ψ).

## Value

Numeric expected log-likelihood value.

## Details

elboTerms - ELBO component functions for VEMSmooth
