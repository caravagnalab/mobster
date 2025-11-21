# Plot summary for model selection.

This plot is usefull to understand the different fits and their rank
with respect to some scoring used to select the best model. The plot
shows alternative solutions and their rank as well.

## Usage

``` r
plot_model_selection(x, TOP = 6, nx = 3, ny = 2, ...)
```

## Arguments

- x:

  A list of fits computed via `mobster_fit`.

- TOP:

  The first `TOP` fits are used (by ranking).

## Value

A complex figure with all plots arranged using both `ggpubr` and
`cowplot`.

## Examples

``` r
data('fit_example', package = 'mobster')
plot_model_selection(fit_example)
#> Required TOP-6 solutions, but only 5 are available.
```
