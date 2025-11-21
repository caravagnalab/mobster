# Summaries for an object of class `'dbpmm'` is like a print.

Summaries for an object of class `'dbpmm'` is like a print.

## Usage

``` r
# S3 method for class 'dbpmm'
print(x, ...)
```

## Arguments

- x:

  An obj of class `'dbpmm'`.

- ...:

## Value

nothing.

## Examples

``` r
data(fit_example)
print(fit_example$best)
#> ── [ MOBSTER ] My MOBSTER model n = 5000 with k = 2 Beta(s) and a tail ─────────
#> ● Clusters: π = 55% [C1], 31% [Tail], and 14% [C2], with π > 0.
#> ● Tail [n = 1370, 31%] with alpha = 1.2.
#> ● Beta C1 [n = 2784, 55%] with mean = 0.48.
#> ● Beta C2 [n = 846, 14%] with mean = 0.15.
#> ℹ Score(s): NLL = -5671.5; ICL = -10359.09 (-11266.35), H = 907.26 (0). Fit
#> converged by MM in 75 steps.
```
