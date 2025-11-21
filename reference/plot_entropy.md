# Plot the entropy of a MOBSTER mixture.

Returns a plot of the entropy and the reduced entropy for this mixture,
binning the domain with bins of size 1e-3. The full entropy is coloured
with a gradient, and the reduced is dashed.

## Usage

``` r
plot_entropy(x)
```

## Arguments

- x:

  An object of class `"dbpmm"`.

## Value

A ggplot object for the plot.

## Examples

``` r
data(fit_example)
plot_entropy(fit_example$best)
```
