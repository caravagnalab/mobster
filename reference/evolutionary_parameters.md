# Extract evolutionary parameters from a MOBSTER fit

Mutation rate, time of emergence and selection coefficient of subclones
are calculated. These values are calculated based on a population
genetics model of tumour evolution see Williams et al. 2016 and 2018 for
more details (Nature Genetics).

The mutation rate scaled by the probability of lineage survival
\\\beta\\, \\\mu/\beta\\ is given by: \$\$\mu/\beta = M / (1/fmin -
1/fmax)\$\$ Where \\fmin\\ is the minimum VAF and \\fmax\\ is the
maximum, and \\M\\ is the number of mutations between \\fmin\\ and
\\fmax\\.

Selection is defined as the relative growth rates of host tumour cell
populations (\\\lambda h\\) vs subclone (\\\lambda s\\): \$\$1+s=\lambda
h / \lambda s\$\$

## Usage

``` r
evolutionary_parameters(
  x,
  Nmax = 10^10,
  lq = 0.1,
  uq = 0.9,
  ploidy = 2,
  ncells = 2
)
```

## Arguments

- x:

  An object fit by MOBSTER.

- Nmax:

  Time when tumour is sampled (in tumour doublings)

- lq:

  lower quantile of VAF (0.05)

- uq:

  upper quantile of VAF (0.95)

- ploidy:

  of mutations (2)

- ncells:

  Number of cells that accumulate mutations at each division 1 or 2,
  default is 1

## Value

Mutation rate, time of emergence and selection coefficient of subclones.

## Examples

``` r
data('fit_example', package = 'mobster')
evolutionary_parameters(fit_example)
#> # A tibble: 1 × 7
#>      mu exponent  time subclonefrequency subclonemutations cluster     s
#>   <dbl>    <dbl> <dbl>             <dbl>             <dbl> <chr>   <dbl>
#> 1  73.5     2.25  5.98             0.298              695. C2      0.177
```
