# Print a Model-Comparison Table

Prints the table
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
return. The default view keeps the columns with the strongest record for
choosing the number of classes, so the table fits one screen: the
log-likelihood and parameter count, `BIC` (the best-supported
information criterion in the simulations of Nylund, Asparouhov & Muthen,
2007), `SABIC` (better than BIC for categorical indicators at smaller
samples, same source), `AIC3` (strong for latent class models with
categorical indicators; Dias, 2006), `Entropy` (classification quality,
not a criterion for K; Masyn, 2013) and the VLMR p-values when they were
requested. `columns = "all"` adds `AIC`, `CAIC`, `AWE`, `ICL`, `BF` and
`cmP`, the full table of Masyn (2013). Nothing is dropped from
`x$fit_table`, which keeps every column at full precision.

## Usage

``` r
# S3 method for class 'mixture_comparison'
print(x, columns = c("default", "all"), ...)
```

## Arguments

- x:

  A `mixture_comparison` object.

- columns:

  `"default"` for the short table, `"all"` for every column.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Details

`AWE` and `cmP` are left out of the default view on purpose. In a large
simulation that included misspecified models, Arch and Nylund-Gibson
(2026) found `BIC` and `CAIC` the most reliable indices, with no
significant difference between the two; `AWE` performed worst overall,
most of all at small samples; and `cmP`, being a rescaling of the `BIC`
differences in the table, tracked `BIC` so closely that printing it adds
no independent evidence and can make one index's choice look like two.
The classes in every simulated population were well separated (entropy
about .84 to .88), so the ranking is firmest there.

A log-likelihood marked `*` is a maximum that only one start reached.
`BF` spans many orders of magnitude, and its reading bands stop at 100
and start at 0.1 (Wagenmakers, 2007), so it is printed as `> 100` and
`< 0.1` beyond them; `cmP` below 0.001 is printed as `< 0.001`.

## References

Arch, D. A. N., & Nylund-Gibson, K. (2026). Unmodeled DIF in latent
class analysis: Consequences for class enumeration. *Methods in
Psychology*, *14*, 100243.
[doi:10.1016/j.metip.2026.100243](https://doi.org/10.1016/j.metip.2026.100243)

Dias, J. G. (2006). Latent class analysis and model selection. In M.
Spiliopoulou et al. (Eds.), *From Data and Information Analysis to
Knowledge Engineering* (pp. 95-102). Springer.

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, *14*(4),
535-569.
[doi:10.1080/10705510701575396](https://doi.org/10.1080/10705510701575396)

Wagenmakers, E.-J. (2007). A practical solution to the pervasive
problems of p values. *Psychonomic Bulletin & Review*, *14*(5), 779-804.
[doi:10.3758/BF03194105](https://doi.org/10.3758/BF03194105)
