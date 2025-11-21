# Filter MOBSTER output clusters.

This function can filter out the clusters computed by MOBSTER based on
two criteria: the mixing proportion value, the number of mutations
assigned and the variance of the Beta clusters.

For all criteria a scalar should be given as input. The return object
will contain only the clusters that pass all filters. If any cluster is
dropped the latent variables are re-computed, as well as the clustering
assignments and the mixing proportions (all mutations will be still
assigned after clusters' removal).

## Usage

``` r
choose_clusters(
  x,
  pi_cutoff = 0.02,
  N_cutoff = 10,
  Beta_variance_cutoff = 1e-04,
  verbose = FALSE
)
```

## Arguments

- x:

  A MOBSTER fit object.

- pi_cutoff:

  The cutoff on the mixing proportions, default is 0.02.

- N_cutoff:

  The cutoff on the number of mutations assigned to a cluster, default
  is 10.

- Beta_variance_cutoff:

  Minimum variance for a Beta peak.

- verbose:

  If outputs should be reported to screen or not, default is no.

## Value

A MOBSTER fit object where clusters are larger than `pi_cutoff` and
contain at least `N_cutoff`. If no such cluster exists an error is
generated.

## Examples

``` r
data('fit_example', package = 'mobster')

# Does not change anything (no filter triggered)
choose_clusters(fit_example$best)
#> ── [ MOBSTER ] My MOBSTER model n = 5000 with k = 2 Beta(s) and a tail ─────────
#> ● Clusters: π = 55% [C1], 31% [Tail], and 14% [C2], with π > 0.
#> ● Tail [n = 1370, 31%] with alpha = 1.2.
#> ● Beta C1 [n = 2784, 55%] with mean = 0.48.
#> ● Beta C2 [n = 846, 14%] with mean = 0.15.
#> ℹ Score(s): NLL = -5671.5; ICL = -10359.09 (-11266.35), H = 907.26 (0). Fit
#> converged by MM in 75 steps.

# Remove one Beta component because it has less than 100 points (renders the fit very poor)
choose_clusters(fit_example$best, N_cutoff = 100)
#> ── [ MOBSTER ] My MOBSTER model n = 5000 with k = 2 Beta(s) and a tail ─────────
#> ● Clusters: π = 55% [C1], 31% [Tail], and 14% [C2], with π > 0.
#> ● Tail [n = 1370, 31%] with alpha = 1.2.
#> ● Beta C1 [n = 2784, 55%] with mean = 0.48.
#> ● Beta C2 [n = 846, 14%] with mean = 0.15.
#> ℹ Score(s): NLL = -5671.5; ICL = -10359.09 (-11266.35), H = 907.26 (0). Fit
#> converged by MM in 75 steps.
```
