# Likelihood-Ratio Test for Two Nested Models (deprecated name)

Deprecated. The test was never specific to longitudinal models — it
accepts any pair of nested fits, and most uses of it in this package are
cross-sectional multiple-group invariance tests. Use
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).

## Usage

``` r
longitudinal_lrt(restricted, full)
```

## Arguments

- restricted:

  The more constrained model (fewer parameters).

- full:

  The less constrained model.

## Value

A list of class `"lr_test"`.

## See also

[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
