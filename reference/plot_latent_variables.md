# Plot the latent variables of the mixture.

It renders a heatmap where the latent variables (reponsibilities) are
shown and colured according to their value. This function also calls
function `Clusters`, using a parameter that determines if a point is not
to be assigned its best cluster based on a cutoff.

## Usage

``` r
plot_latent_variables(x, cutoff_assignment = 0)
```

## Arguments

- x:

  A MOBSTER fit.

- cutoff_assignment:

  The parameter used to call function `Clusters`, which does not assign
  a point to its best cluster if the value of the corresponding latent
  variable is not above the cutoff.

## Value

A plot of the latent variables of the mixture.

## Examples

``` r
data('fit_example', package = 'mobster')
plot_latent_variables(fit_example$best)

plot_latent_variables(fit_example$best, cutoff_assignment = .9)
```
