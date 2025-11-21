# Return a tabular representation of a model.

This functionr returns a table with all the parameters fits (one per
column), the scores of the model, and the SSE of the fit versus data.

## Usage

``` r
to_string(x)
```

## Arguments

- x:

  A MOBSTER fit.

## Value

## Examples

``` r
data('fit_example', package = 'mobster')

to_string(fit_example$best)
#> # A tibble: 1 × 39
#>   K_beta Mean_C1 Mean_C2 Mean_Tail     N  N_C1  N_C2 N_Tail Scale_Tail
#> *  <int>   <dbl>   <dbl>     <dbl> <int> <dbl> <dbl>  <dbl>      <dbl>
#> 1      2   0.478   0.149     0.251  5000  2784   846   1370     0.0500
#> # ℹ 30 more variables: Shape_Tail <dbl>, Variance_C1 <dbl>, Variance_C2 <dbl>,
#> #   Variance_Tail <dbl>, pi_C1 <dbl>, pi_C2 <dbl>, pi_Tail <dbl>, rcc_C1 <lgl>,
#> #   rcc_C2 <lgl>, rcc_Tail <lgl>, tail <lgl>, NLL <dbl>, BIC <dbl>, AIC <dbl>,
#> #   entropy <dbl>, ICL <dbl>, reduced.entropy <dbl>, reICL <dbl>, size <dbl>,
#> #   sse_total <dbl>, sse_0_0.1 <dbl>, sse_0.1_0.2 <dbl>, sse_0.2_0.3 <dbl>,
#> #   sse_0.3_0.4 <dbl>, sse_0.4_0.5 <dbl>, sse_0.5_0.6 <dbl>, sse_0.6_0.7 <dbl>,
#> #   sse_0.7_0.8 <dbl>, sse_0.8_0.9 <dbl>, sse_0.9_1 <dbl>
```
