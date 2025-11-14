# devPsi - Derivatives of Correlation Matrix with respect to w

Computes the first and second derivatives of the Ψ correlation matrix
with respect to the decay parameter w.

## Usage

``` r
devPsi(Xt, Xp, w = 1)
```

## Arguments

- Xt:

  Numeric vector of time points.

- Xp:

  Numeric vector of time points.

- w:

  Numeric decay parameter.

## Value

A list with:

- dPsi:

  First derivative matrix (∂Ψ/∂w)

- d2Psi:

  Second derivative matrix (∂²Ψ/∂w²)
