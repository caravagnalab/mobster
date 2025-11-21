# Plot the initial density of a fit.

Plot the initial density of a fit.

## Usage

``` r
plot_init(x, colors = c(Tail = "gainsboro"))
```

## Arguments

- x:

  An object of class `"dbpmm"`.

- colors:

  If provided, these colours will be used for each cluster. If a subset
  of colours is provided, palette Set1 from `RColorBrewer` is used. By
  default the tail colour is provided as 'gainsboro'.

## Value

A ggplot object for the plot.

## Examples

``` r
data(fit_example)
plot_init(fit_example$best)
```
