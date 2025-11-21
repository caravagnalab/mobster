# Run a dN/dS analysis on MOBSTER clusters.

This function takes a MOBSTER fit and runs \`dndscv\`
(https://github.com/im3sanger/dndscv) to calculate dN/dS values per
cluster. It computes global dN/dS and per gene dN/dS values and makes a
plot. dN/dS values are computed with the best fitting MOBSTER model.

## Usage

``` r
dnds(
  x,
  mapping = NULL,
  gene_list = NULL,
  colors = c(Tail = "gray"),
  refdb = "hg19",
  dndscv_plot = c("wall", "wmis", "wnon", "wspl", "wtru"),
  ...
)
```

## Arguments

- x:

  A MOBSTER fit object.

- mapping:

  The groups used to compute this statistics are defined by this
  variable. If \`mapping = c(\`A\` = 'G1', \`B\` = 'G1', \`C\` =
  'G2')\`, then mutations from clusters \`A\` and \`B\` will be pooled
  into one group (\`G1\`), while mutations from cluster \`C\` will
  constitute a group themselves. By default, with \`mapping = NULL\`,
  each cluster is a group.

- gene_list:

  An optional vector of gene names to infer dN/dS values, default
  (\`NULL\`) is to use `dndscv` default (whole-exome. This package
  provides lists genes that can be used for this value (essential genes,
  cancer genes, etc.); see package data.

- colors:

  If provided, these colours will be used for each cluster. If a subset
  of colours is provided, palette Set1 from `RColorBrewer` is used. By
  default the tail colour is provided as 'gainsboro'.

- refdb:

  The genome referene to use, default is to use hg19. Other references
  are available from https://github.com/im3sanger/dndscv_data

- dndscv_plot:

  What of the dndscv scores should be visualized in a plot, by default
  all the statistcs are reported. One can use \`dndscv_plot = wall\` to
  get only the global dnds value.

- ...:

  Extra parameters forwarded to `dndscv`.

## Value

The fit object is a list with the summary table and the observation
counts reported by package `dndscv`, together with a `ggplot` plot for
the results.

## Examples

``` r
# Example run with real data
data('LUFF76_lung_sample', package = 'mobster')

clusters = Clusters(LUFF76_lung_sample$best)

dnds_stats = dnds(clusters, gene_list = NULL)
#> Missing 'sample' column, assuming mutations from a single patient (adding a sample label otherwise).
#> Warning: Unknown or uninitialised column: `sample`.
#> ℹ 2298 mutations; 'by cluster' groups in 0 samples, with no genes (default dndscv).
#> [refdb = hg19] Removing chr from chromosome names for hg19 reference compatability
#> 
#>   C1 Tail 
#> 1222 1076 
#> 
#> ── Running dndscv ──────────────────────────────────────────────────────────────
#> 
#> ── Group Tail 
#> [1] Loading the environment...
#> [2] Annotating the mutations...
#> Warning: Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.
#> Loading required namespace: GenomeInfoDb
#> [3] Estimating global rates...
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> [4] Running dNdSloc...
#> [5] Running dNdScv...
#>     Regression model for substitutions (theta = 6.69e-05).
#> 
#> ── Group C1 
#> [1] Loading the environment...
#> [2] Annotating the mutations...
#> Warning: Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.
#> [3] Estimating global rates...
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted rates numerically 0 occurred
#> [4] Running dNdSloc...
#> [5] Running dNdScv...
#> dndscv error
#> Error in while ((it <- it + 1) < limit && abs(del) > eps) {: missing value where TRUE/FALSE needed
#> 
#> ── dndscv results ────────────────────────────── wall, wmis, wnon, wspl, wtru ──
#> # A tibble: 5 × 5
#>   name            mle cilow cihigh dnds_group
#>   <chr>         <dbl> <dbl>  <dbl> <chr>     
#> 1 wmis  0.731         0.113   4.73 Tail      
#> 2 wnon  0.00000000804 0     Inf    Tail      
#> 3 wspl  0.0000000260  0     Inf    Tail      
#> 4 wtru  0.0000000170  0     Inf    Tail      
#> 5 wall  0.701         0.109   4.51 Tail      
```
