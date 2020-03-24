
# mobster <a href='https://caravagn.github.io/mobster'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/caravagn/mobster.svg?branch=master)](https://travis-ci.org/caravagn/mobster)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/badge/Part%20of-evoverse-blue.svg)](https://caravagn.github.io/evoverse)
<!-- badges: end -->

`mobster` is a package that implements a model-based clustering approach
to subclonal deconvolution from cancer genome sequencing data
([Caravagna et al.;
https://doi.org/10.1101/586560](https://www.biorxiv.org/content/10.1101/586560v1)).

The paper as of today is still under review.

The package integrates evolutionary theory and Machine-Learning to
analyze bulk sequencing data of a cancer biopsy - ideally,
high-resolution whole-genome sequencing data (e.g., WGS \>100x).
`mobster` fits can be computed via moment-matching or
maximum-likelihood, and model selection can be done with multiple
likelihood-based scores (BIC, ICL and reICL, a new reduced-entropy
variation to ICL). S3 objects allow easy visualization of the data, the
fits and the quality of the model. Parametric and nonparametric
bootstrap routines are available to assess confidence of the parameters.

`mobster` is part of the `evoverse` set of [R
packages](https://caravagn.github.io/evoverse) to implement Cancer
Evolution
analyses.

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagn.github.io/mobster/-steelblue.svg)](https://caravagn.github.io/mobster)

### Installation

You can install the released version of `mobster` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagn/mobster")
```

-----

#### Copyright and contacts

Giulio Caravagna, PhD. *Institute of Cancer Research, London,
UK*.

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-informational.svg?style=social)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/caravagn-informational.svg?style=social&logo=GitHub)](https://github.com/caravagn)
[![](https://img.shields.io/badge/@gcaravagna-informational.svg?style=social&logo=Twitter)](https://twitter.com/gcaravagna)
[![](https://img.shields.io/badge/Homepage-informational.svg?style=social&logo=Google)](https://sites.google.com/site/giuliocaravagna/)
