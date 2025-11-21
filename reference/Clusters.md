# Return the data with the hard clustering assigments.

Extract the clustering assignments, and return the data. Assignments can
be computed subsetting the mutations that have posterior values
(reponsibilities) below a certain cutoff (default `0` - all assigments);
non-assigned mutations have `NA` as cluster label.

## Usage

``` r
Clusters(x, cutoff_assignment = 0)
```

## Arguments

- x:

  A MOBSTER fit.

- cutoff_assignment:

  Cutoff to compute hard clustering assignments.

## Value

The data stored in `x$data` with a column `label` reporting the assigned
cluster, or `NA` if the maximum cluster probability is below the
threshold value `cutoff_assignment`.

## Examples

``` r
data('fit_example', package = 'mobster')
Clusters(fit_example$best)
#> # A tibble: 5,000 × 5
#>      VAF cluster    Tail    C1       C2
#>    <dbl> <chr>     <dbl> <dbl>    <dbl>
#>  1 0.497 C1      0.00736 0.993 5.22e-27
#>  2 0.490 C1      0.00669 0.993 4.42e-26
#>  3 0.470 C1      0.00705 0.993 1.31e-23
#>  4 0.517 C1      0.0130  0.987 1.83e-29
#>  5 0.506 C1      0.00903 0.991 3.86e-28
#>  6 0.440 C1      0.0179  0.982 9.68e-20
#>  7 0.428 C1      0.0347  0.965 3.88e-18
#>  8 0.523 C1      0.0164  0.984 3.97e-30
#>  9 0.482 C1      0.00648 0.994 3.87e-25
#> 10 0.499 C1      0.00759 0.992 3.20e-27
#> # ℹ 4,990 more rows

# Add some cutoff to filter assignments
Clusters(fit_example$best, cutoff_assignment = .8)
#> # A tibble: 5,000 × 5
#>      VAF cluster    Tail    C1       C2
#>    <dbl> <chr>     <dbl> <dbl>    <dbl>
#>  1 0.497 C1      0.00736 0.993 5.22e-27
#>  2 0.490 C1      0.00669 0.993 4.42e-26
#>  3 0.470 C1      0.00705 0.993 1.31e-23
#>  4 0.517 C1      0.0130  0.987 1.83e-29
#>  5 0.506 C1      0.00903 0.991 3.86e-28
#>  6 0.440 C1      0.0179  0.982 9.68e-20
#>  7 0.428 C1      0.0347  0.965 3.88e-18
#>  8 0.523 C1      0.0164  0.984 3.97e-30
#>  9 0.482 C1      0.00648 0.994 3.87e-25
#> 10 0.499 C1      0.00759 0.992 3.20e-27
#> # ℹ 4,990 more rows
```
