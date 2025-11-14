# checkConvergence - ELBO Convergence test

Checks whether the algorithm has converged based on ELBO changes.

## Usage

``` r
checkConvergence(elbo_c, elbo_prev, convergence_threshold)
```

## Arguments

- elbo_c:

  Current ELBO value.

- elbo_prev:

  Previous ELBO value.

- convergence_threshold:

  Numeric convergence threshold.

## Value

Logical; TRUE if converged.
