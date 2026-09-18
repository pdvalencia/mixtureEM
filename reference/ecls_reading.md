# Reading proficiency from kindergarten to first grade (ECLS-K)

Five binary reading-proficiency indicators measured on the same 3,575
children at four occasions – the fall and spring of kindergarten and the
fall and spring of first grade – from the Early Childhood Longitudinal
Study, Kindergarten Class of 1998-99 (ECLS-K), together with a household
poverty indicator. The five indicators mark mastery of successive
reading skills, which is what makes the data a textbook case for a
stage-sequential latent transition analysis (Kaplan, 2008): the latent
statuses are stages of learning to read, and children are expected to
move forward through them. Used in
[`vignette("lta")`](https://pdvalencia.github.io/mixtureEM/articles/lta.md)
and
[`vignette("rilta")`](https://pdvalencia.github.io/mixtureEM/articles/rilta.md).

## Usage

``` r
ecls_reading
```

## Format

A data frame with 3,575 rows (children) and 21 columns. The indicator
columns are in the wide, occasion-major layout
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
reads by default: the five items at `t1`, then the five at `t2`, and so
on.

- poverty:

  Household below the U.S. census poverty threshold (0/1; 19% are 1).

- letters_t1, letters_t2, letters_t3, letters_t4:

  Letter recognition mastered (0/1).

- beginning_t1, ..., beginning_t4:

  Beginning sounds mastered (0/1).

- ending_t1, ..., ending_t4:

  Ending sounds mastered (0/1).

- sight_t1, ..., sight_t4:

  Sight words mastered (0/1).

- context_t1, ..., context_t4:

  Words in context mastered (0/1).

Occasions: `t1` = fall of kindergarten, `t2` = spring of kindergarten,
`t3` = fall of first grade, `t4` = spring of first grade. There are no
missing values.

## Source

U.S. Department of Education, National Center for Education Statistics.
Early Childhood Longitudinal Study, Kindergarten Class of 1998-99,
kindergarten-first grade public-use child file.
<https://nces.ed.gov/ecls/dataproducts.asp>. ECLS-K public-use data are
produced by a U.S. federal agency and are in the public domain (17
U.S.C. sec. 105); no permission is required to use or redistribute them.
The analytic extract shipped here defines the indicators as in Kaplan
(2008); see `data-raw/ecls_reading.R`.

## References

Kaplan, D. (2008). An overview of Markov chain methods for the study of
stage-sequential developmental processes. *Developmental Psychology*,
*44*(2), 457-467.
[doi:10.1037/0012-1649.44.2.457](https://doi.org/10.1037/0012-1649.44.2.457)

## Examples

``` r
data(ecls_reading)
# Mastery rises from occasion to occasion, and in the order of the skills
round(colMeans(ecls_reading[, -1]), 2)
#>   letters_t1 beginning_t1    ending_t1     sight_t1   context_t1   letters_t2 
#>         0.64         0.32         0.18         0.03         0.01         0.90 
#> beginning_t2    ending_t2     sight_t2   context_t2   letters_t3 beginning_t3 
#>         0.73         0.53         0.15         0.05         0.93         0.82 
#>    ending_t3     sight_t3   context_t3   letters_t4 beginning_t4    ending_t4 
#>         0.68         0.27         0.10         0.97         0.93         0.90 
#>     sight_t4   context_t4 
#>         0.81         0.45 
```
