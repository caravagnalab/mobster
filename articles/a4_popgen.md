# 4. Population Genetics statistics

``` r
library(mobster)
library(tidyr)
library(dplyr)
```

Population Genetics statistics can be extracted from a MOBSTER model.

``` r
data('fit_example', package = 'mobster')
print(fit_example$best)
#> ── [ MOBSTER ] My MOBSTER model n = 5000 with k = 2 Beta(s) and a tail ─────────
#> ● Clusters: π = 55% [C1], 31% [Tail], and 14% [C2], with π > 0.
#> ● Tail [n = 1370, 31%] with alpha = 1.2.
#> ● Beta C1 [n = 2784, 55%] with mean = 0.48.
#> ● Beta C2 [n = 846, 14%] with mean = 0.15.
#> ℹ Score(s): NLL = -5671.5; ICL = -10359.09 (-11266.35), H = 907.26 (0). Fit
#> converged by MM in 75 steps.

evolutionary_parameters(fit_example)
#> # A tibble: 1 × 7
#>      mu exponent  time subclonefrequency subclonemutations cluster     s
#>   <dbl>    <dbl> <dbl>             <dbl>             <dbl> <chr>   <dbl>
#> 1  73.5     2.25  5.98             0.298              695. C2      0.177
```

The mutation rate `mu` (cell division units) scaled by the probability
of lineage survival $\beta$, $\mu/\beta$, is given by:
$$\mu/\beta = \frac{M}{\left( \frac{1}{f_{\text{min}}} - \frac{1}{f_{\text{max}}} \right)}$$
Where $f_{\text{min}}$ is the minimum VAF and $f_{\text{max}}$ is the
maximum, and $M$ is the number of mutations between $f_{\text{min}}$ and
$f_{\text{max}}$.

Selection is defined as the relative growth rates of host tumour cell
populations ($\lambda h$) vs subclone ($\lambda s$):
$$1 + s = \frac{\lambda h}{\lambda s}$$

The mathematical details of these computations are described in the main
paper, and baesd on the population genetics model of tumour evolutionin
Williams et al. 2016 and 2018 (Nature Genetics).
