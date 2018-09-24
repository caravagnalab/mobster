# MOBSTER

Model-based subclonal deconvolution from bulk sequencing. The statistical model uses a finite Dirichlet mixture model with Beta and Pareto components. `K` Beta random variables model `K` subclones in the VAF distribution, and one Pareto power law tail predicted by Population Genetics, describes alleles under neutral evolution. Fits are via moment-matching or MLE.
