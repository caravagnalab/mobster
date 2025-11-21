# Package index

## Available data

Data released with the package

- [`fit_example`](https://caravagnalab.github.io/mobster/reference/fit_example.md)
  : Example MOBSTER fit.
- [`LU4_lung_sample`](https://caravagnalab.github.io/mobster/reference/LU4_lung_sample.md)
  : MOBSTER fit for the LU4 lung sample.
- [`LUFF76_lung_sample`](https://caravagnalab.github.io/mobster/reference/LUFF76_lung_sample.md)
  : MOBSTER fit for the LUFF76 lung sample.
- [`cancer_genes_dnds`](https://caravagnalab.github.io/mobster/reference/cancer_genes_dnds.md)
  : List of cancer genes to compute dnds values.

## Fitting functions

Functions to cluster the data; density of a `mobster` fit and sampling
functions.

- [`mobster_fit()`](https://caravagnalab.github.io/mobster/reference/mobster_fit.md)
  : Fit a model with MOBSTER.
- [`choose_clusters()`](https://caravagnalab.github.io/mobster/reference/choose_clusters.md)
  : Filter MOBSTER output clusters.
- [`Clusters()`](https://caravagnalab.github.io/mobster/reference/Clusters.md)
  : Return the data with the hard clustering assigments.
- [`Clusters_denovo()`](https://caravagnalab.github.io/mobster/reference/Clusters_denovo.md)
  : Assign new observations to the clusters inside a MOBSTER fit.
- [`ddbpmm()`](https://caravagnalab.github.io/mobster/reference/ddbpmm.md)
  : Density function for MOBSTER.
- [`rdbpmm()`](https://caravagnalab.github.io/mobster/reference/rdbpmm.md)
  : Generate a random sample from a MOBSTER model.
- [`random_dataset()`](https://caravagnalab.github.io/mobster/reference/random_dataset.md)
  : Generate a random MOBSTER model and data.
- [`to_string()`](https://caravagnalab.github.io/mobster/reference/to_string.md)
  : Return a tabular representation of a model.

## S3 functions

The `mobster` S3 object.

- [`print(`*`<dbpmm>`*`)`](https://caravagnalab.github.io/mobster/reference/print.dbpmm.md)
  :

  Summaries for an object of class `'dbpmm'` is like a print.

- [`summary.dbpmm()`](https://caravagnalab.github.io/mobster/reference/summary.dbpmm.md)
  :

  Summary for an object of class `'dbpmm'` is a print.

- [`plot.dbpmm()`](https://caravagnalab.github.io/mobster/reference/plot.dbpmm.md)
  : Plot a MOBSTER fit.

## Plotting functions

These functions can be used to plot data and `mobster` fits.

- [`plot_NLL()`](https://caravagnalab.github.io/mobster/reference/plot_NLL.md)
  : Plot the NLL trace.
- [`plot_entropy()`](https://caravagnalab.github.io/mobster/reference/plot_entropy.md)
  : Plot the entropy of a MOBSTER mixture.
- [`plot_fit_scores()`](https://caravagnalab.github.io/mobster/reference/plot_fit_scores.md)
  : Plot the scores for model selection.
- [`plot_gofit()`](https://caravagnalab.github.io/mobster/reference/plot_gofit.md)
  : Plot the goodness of fit.
- [`plot_init()`](https://caravagnalab.github.io/mobster/reference/plot_init.md)
  : Plot the initial density of a fit.
- [`plot_latent_variables()`](https://caravagnalab.github.io/mobster/reference/plot_latent_variables.md)
  : Plot the latent variables of the mixture.
- [`plot_mixing_proportions()`](https://caravagnalab.github.io/mobster/reference/plot_mixing_proportions.md)
  : Plot the mixing proportions of the mixture.
- [`plot_model_selection()`](https://caravagnalab.github.io/mobster/reference/plot_model_selection.md)
  : Plot summary for model selection.

## Bootstrap functions

Bootstrap functions for `mobster` fits, and plotting functions.

- [`mobster_bootstrap()`](https://caravagnalab.github.io/mobster/reference/mobster_bootstrap.md)
  : Bootstrap a MOBSTER fit.
- [`bootstrapped_statistics()`](https://caravagnalab.github.io/mobster/reference/bootstrapped_statistics.md)
  : Compute boostrap statistics from a bootstrap run.
- [`plot_bootstrap_Beta()`](https://caravagnalab.github.io/mobster/reference/plot_bootstrap_Beta.md)
  : Plot the boostrapped tail parameters.
- [`plot_bootstrap_coclustering()`](https://caravagnalab.github.io/mobster/reference/plot_bootstrap_coclustering.md)
  : Plot the boostrapped co-clustering probability.
- [`plot_bootstrap_mixing_proportions()`](https://caravagnalab.github.io/mobster/reference/plot_bootstrap_mixing_proportions.md)
  : Plot the boostrapped mixing proportions
- [`plot_bootstrap_model_frequency()`](https://caravagnalab.github.io/mobster/reference/plot_bootstrap_model_frequency.md)
  : Plot the boostrapped model frequency
- [`plot_bootstrap_tail()`](https://caravagnalab.github.io/mobster/reference/plot_bootstrap_tail.md)
  : Plot the boostrapped tail parameters.

## Extra features

Other analyses for `mobster` fits.

- [`dnds()`](https://caravagnalab.github.io/mobster/reference/dnds.md) :
  Run a dN/dS analysis on MOBSTER clusters.
- [`evolutionary_parameters()`](https://caravagnalab.github.io/mobster/reference/evolutionary_parameters.md)
  : Extract evolutionary parameters from a MOBSTER fit
- [`get_clone_trees()`](https://caravagnalab.github.io/mobster/reference/get_clone_trees.md)
  : Return clone trees from the fit.
- [`load_vcf()`](https://caravagnalab.github.io/mobster/reference/load_vcf.md)
  : Load data from a VCF file
