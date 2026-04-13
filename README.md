# Conditional Mean Independence Simulations

This repository contains three R scripts implementing simulation studies for a spline-based conditional mean independence test with wild bootstrap calibration.

The examples share a common workflow:

1. Generate data under a null or alternative model.
2. Build a spline basis for the conditioning variable `Z`.
3. Regress `Y` on the basis expansion of `Z`.
4. Compute the MDD-based test statistic using the residuals and `(Z, X)`.
5. Use wild bootstrap to estimate critical values and p-values.
6. Repeat over many Monte Carlo replications to estimate empirical size and power.

## Files

- `CMI_Example 5.1.R`
  A linear mean function for `E(Y | Z)` with
  `Z ~ Unif(-1,1)`, `X ~ Unif(-1,1)`, and
  `Y = Z + 0.1 * eps + r * X`.

- `CMI_Example 5.2.R`
  A nonlinear mean function for `E(Y | Z)`,
  while Z and X are linearly related with
  `Y = (Z^2 - 1) + 0.1 * eps + r * X`.

- `CMI_Example 5.3.R`
  An example adopted from Park et al (2015) where `Y` and `X` are dependent with
  `X = Z + U2` and
  `Y = Z + U1 + r * X`.

## Main Components

Each script includes the following helper functions:

- `matpower()`
  Computes matrix powers used to standardize the spline basis.

- `fit_with_qr()`
  Uses QR decomposition for least-squares fitting without repeated `lm()` overhead inside the bootstrap loop.

- `generate_sample()`
  Defines the data-generating mechanism for the corresponding example.

- `build_basis()`
  Builds either a B-spline basis (`bs`) or natural spline basis (`ns`) for `Z`, then centers and standardizes it.

- `select_p()`
  Optionally selects the number of basis functions using:
  - generalized cross-validation (`p_method = "gcv"`), or
  - a Mallows-type information criterion (`p_method = "mallows"`).

- `run_case()`
  Runs one Monte Carlo replication and returns:
  - rejection indicators at levels `0.10`, `0.05`, `0.01`
  - bootstrap p-value
  - test statistic
  - bootstrap critical values
  - selected basis size

## Dependencies

Install the required R packages before running the scripts:

```r
install.packages(c("EDMeasure", "splines", "MASS"))
```

Notes:

- `splines` is part of base R, but loading it explicitly is still helpful.
- `MASS` is only needed for `CMI_Example 5.2.R`.

## Example Usage

Run a single replication with a fixed number of basis functions:

```r
run_case(
  p = 10,
  i = 1,
  n = 100,
  p_method = "fixed",
  basis_type = "bs",
  degree = 3,
  boundary_knots = c(-1, 1)
)
```

Run a single replication with automatic basis selection by GCV:

```r
run_case(
  i = 1,
  n = 100,
  p_method = "gcv",
  p_candidates = 3:15,
  basis_type = "bs",
  degree = 3,
  boundary_knots = c(-1, 1)
)
```

## Simulation Output

Each script runs 500 Monte Carlo replications for selected values of `p` and stores:

- `results_by_p`
  A list of `500 x 3` matrices of rejection indicators.

- `sizes`
  Empirical rejection rates under the null.

- `powers`
  Empirical rejection rates under the alternative.

Rows correspond to different basis sizes, and columns correspond to significance levels:

- `alpha_0.10`
- `alpha_0.05`
- `alpha_0.01`

## Remarks

- The basis size can be fixed manually or selected automatically.
- The scripts currently use wild bootstrap with `BB = 499`.
- Boundary knots differ by example and should be chosen to match the support of `Z`.

## Reproducibility

- Each replication uses `set.seed(i)` inside `generate_sample()`.
- Bootstrap draws use deterministic seeds of the form `set.seed(100 + k)`.
- Results are therefore reproducible for a given script, seed, and parameter setting.
