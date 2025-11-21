# Plot the goodness of fit.

Plot the SSE (sum of squared error) as a proxy for the goodness of fit.
For the `TOP` available solutions the SSE trace is shown.

## Usage

``` r
plot_gofit(x, TOP = 6, nx = 3, ny = 2)
```

## Arguments

- x:

  A list of fits computed via `mobster_fit`.

- TOP:

  The first `TOP` fits are used.

- nx:

  Columns in the matrix layout.

- ny:

  Rows in the matrix layout.

## Value

A plot of the goodness of fit.

## Examples

``` r
data('fit_example', package = 'mobster')
plot_gofit(fit_example, TOP = 3)
```
