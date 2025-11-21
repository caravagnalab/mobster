# Plot the scores for model selection.

Plots the scores via ICL, reICL, BIC and AIC which can be used for model
selection. It allows to easily see if the model selected as best is
consistently better for all scores.

## Usage

``` r
plot_fit_scores(x)
```

## Arguments

- x:

  A list of fits computed via `mobster_fit`.

## Value

A ggplot figure with the scores for model selection.

## Examples

``` r
data('fit_example', package = 'mobster')
plot_fit_scores(fit_example)
```
