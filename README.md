
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epiGWAS

[![Rdoc](http://www.rdocumentation.org/badges/version/epiGWAS)](http://www.rdocumentation.org/packages/epiGWAS)
[![Travis build
status](https://travis-ci.org/EpiSlim/epiGWAS.svg?branch=master)](https://travis-ci.org/EpiSlim/epiGWAS)

This package implements a number of methods for detecting pure epistatic
interactions with a predetermined target variant. The common denominator
lies in the use of propensity scores to filter out the main effects of
the rest of genotype. The methods incorporate propensity scores in two
different ways: either in the sample weights (outcome weighted learning)
or in the response (modified outcome).

## Installation

You can install the released version of epiGWAS from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("epiGWAS")
```

The latest development version is directly available from
[GitHub](https://github.com):

``` r
install.packages("devtools")
devtools::install_github("EpiSlim/epiGWAS")
```

## Usage Examples

The example below illustrates how to use our methods on a synthetic
dataset:

``` r
require("epiGWAS")

# Genotype simulation
set.seed(542)
n_samples <- 300
p <- 450
# Genotypes matrix with {0, 1, 2} SNP values
genotypes <- matrix(
  (runif(n_samples * p, min = 0, max = 1) <
     runif(n_samples * p, min = 0, max = 1)) +
    (runif(n_samples * p, min = 0, max = 1) <
       runif(n_samples * p, min = 0, max = 1)),
  ncol = p, nrow = n_samples, dimnames = list(NULL, paste0("SNP_", seq_len(p)))
)

# Phenotype simulation
target <- "SNP_56"
syner <- paste0("SNP_", sample.int(p, 10))
size_effects <- rnorm(10) 
binarized <- genotypes[, target] > 1
risk <-   (2 * binarized - 1) * (genotypes[, syner] %*% size_effects)
risk <- risk - mean(risk) # Centering to balance cases and controls
phenotype <- runif(n_samples) < 1/(1+exp(-risk)) # Logistic model
```

The propensity scores can be estimated using the fastPHASE hidden Markov
model. Make sure to download the
[fastPHASE](http://scheet.org/software.html) executable before running
the `fast_HMM` function.

``` r
hmm <- fast_HMM(genotypes, fp_path = "/path/to/fastPHASE",
  n_state = 4, n_iter = 10)
propensity <- cond_prob(genotypes, target, hmm, binary = FALSE)
propensity <- propensity[cbind(seq(dim(genotypes)[1]), binarized + 1)]
```

All the pieces are now in place to apply our epistasis detection methods
via the `epiGWAS`
function.

``` r
stability_scores <- epiGWAS(binarized, genotypes[, colnames(genotypes) != target], phenotype,
                            propensity, methods = c("OWL", "modified_outcome", "shifted_outcome",
                                        "normalized_outcome", "robust_outcome"), parallel = FALSE)
```

## References

Slim, L., Chatelain, C., Azencott, C.-A., & Vert, J.-P. (2018). Novel
Methods for Epistasis Detection in Genome-Wide Association Studies.
BioRxiv. Retrieved from
<http://biorxiv.org/content/early/2018/10/14/442749>
