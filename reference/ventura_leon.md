# Infidelity behavior patterns in Peruvian young adults

Sixteen binary indicators from the Multidimensional Infidelity
Inventory-S (MII-S; Romero-Palencia et al. 2007) plus six
sociodemographic variables, from a latent class analysis of
infidelity-related thoughts and behaviors in a sample of 400 Peruvian
young adults. Reproduces the four-class solution (infidelity, sexual
desire, affective interest, fidelity) of Ventura-Leon et al. (2025). The
original paper related classes to covariates by assigning each
respondent to their most likely class and then running chi-square/ANOVA
tests on the observed class labels – a classify-and-analyze approach
that ignores classification uncertainty.
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)'s
`predictors=` argument with `n_steps = 3` and ML bias correction is used
instead in the package vignette, which accounts for that uncertainty
when testing the same covariates (sex, age, sexual orientation,
relationship duration).

## Usage

``` r
ventura_leon
```

## Format

A data frame with 400 rows (respondents) and 22 columns:

- sex:

  Factor: `"Male"`, `"Female"`.

- age:

  Age in years.

- birthplace:

  Factor: `"Lima"`, `"Outside Lima"`.

- sexual_orientation:

  Factor: `"Heterosexual"`, `"Not heterosexual"`.

- relationship_duration:

  Factor: `"Short"`, `"Long"`, dichotomized at the sample median (50
  months), per the original study.

- relationship_type:

  Factor: `"Married"`, `"Cohabiting"`, `"In love"`, `"Engaged"`,
  `"Dating"`.

- flirting:

  "I have flirted with other people besides my partner." (0/1)

- romantic_partners:

  "I have had another romantic partner(s)." (0/1)

- emotional_bond:

  "I have formed an emotional bond with someone else besides my
  partner." (0/1)

- romantic_involvement:

  "I have been romantically involved with another person/other people."
  (0/1)

- loved_another:

  "I have loved another person(s) besides my partner." (0/1)

- in_love:

  "I have fallen in love with someone else besides my partner." (0/1)

- thoughts:

  "I have thought about someone else besides my partner." (0/1)

- interest:

  "I have been interested in another person(s) besides my partner."
  (0/1)

- sexual_relations:

  "I have had sexual relations with another person/other people besides
  my partner." (0/1)

- sexual_contact:

  "I have had sexual contact with someone else besides my partner."
  (0/1)

- desired_relations:

  "I have desired to have sexual relations with another person(s)
  besides my partner." (0/1)

- desired_contact:

  "I have desired to have sexual contact with another person(s) besides
  my partner." (0/1)

- sexual_fantasies:

  "I have wanted to fulfill my sexual fantasies with someone else
  besides my partner." (0/1)

- attraction:

  "I have felt attracted to another person(s) besides my partner." (0/1)

- had_sex:

  "I have had sex with another person/other people besides my partner."
  (0/1)

- desired_sex:

  "I have desired to have sex with another person(s) besides my
  partner." (0/1)

## Source

OSF repository <https://osf.io/8csr9/>. The OSF project has no license
flag set; Pablo D. Valencia, a co-author of the source paper and of this
package, has authorized bundling this dataset in mixtureEM on behalf of
the author team. See `data-raw/ventura_leon.R` for the exact recoding
script.

## References

Ventura-Leon, J., Reyes, A., Valencia, P. D., Tocto-Munoz, S.,
Gamboa-Melgar, G., Ruiz-Castro, J., & Lino-Cruz, C. (2025). Exploring
infidelity behavior patterns in a sample of Peruvian young adults: A
latent class analysis. *Journal of Marital and Family Therapy*, *51*,
e70066. [doi:10.1111/jmft.70066](https://doi.org/10.1111/jmft.70066) .
Tables 1-4 report the sociodemographic, item-response and covariate
numbers this recoding and the package vignette are checked against.
