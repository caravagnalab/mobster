# Plot the NLL trace.

It plots the Negative log-Likelihood (NLL) trace of a fit.

## Usage

``` r
plot_NLL(x)
```

## Arguments

- x:

  A MOBSTER fit.

## Value

A plot of the the NLL trace for the input fit.

## Examples

``` r
data('fit_example', package = 'mobster')
plot_NLL(fit_example$best)
```
