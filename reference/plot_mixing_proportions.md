# Plot the mixing proportions of the mixture.

Plot the mixing proportions of the mixture.

## Usage

``` r
plot_mixing_proportions(x, colors = c(Tail = "gainsboro"))
```

## Arguments

- x:

  An object of class 'dbpmm'.

- colors:

  If provided, these colours will be used for each cluster. If a subset
  of colours is provided, palette Set1 from `RColorBrewer` is used. By
  default the tail colour is provided as 'gainsboro'.

## Value

A plot of the mixing proportions of the mixture.

## Examples

``` r
data(fit_example)
plot_mixing_proportions(fit_example$best)
#> Warning: All aesthetics have length 1, but the data has 3 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
```
