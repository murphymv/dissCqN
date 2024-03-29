
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dissCqN

<!-- badges: start -->

[![Repo
Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg?label=Lifecycle)](https://lifecycle.r-lib.org/articles/stages.html)
[![Licence](https://img.shields.io/badge/License-GPL3-green.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![GitHub language
count](https://img.shields.io/github/languages/count/murphymv/dissCqN?label=Languages)
[![R-CMD-check](https://github.com/murphymv/dissCqN/workflows/R-CMD-check/badge.svg)](https://github.com/murphymv/dissCqN/actions)

[![CRAN](https://www.r-pkg.org/badges/version/dissCqN?color=blue)](https://CRAN.R-project.org/package=dissCqN)
![Downloads:
Total](http://cranlogs.r-pkg.org/badges/grand-total/dissCqN)
![Downloads: Monthly](https://cranlogs.r-pkg.org/badges/dissCqN)

[![Donate](https://img.shields.io/badge/PayPal-Donate%20to%20Author-yellow.svg)](https://paypal.me/murphymv1)

<!-- badges: end -->

`dissCqN` is a small package designed to make the process of calculating
multiple or pairwise assemblage dissimilarity via the C<sub>*qN*</sub>
generalisation of similarity indices (Chao et al., 2008; Jost et al.,
2011) relatively straightforward and fast. Although C<sub>*qN*</sub> can
also be calculated using the `SpadeR` package
(e.g. [`SpadeR::SimilarityMult()`](https://rdrr.io/cran/SpadeR/man/SimilarityMult.html)
and
[`SpadeR::SimilarityPair()`](https://rdrr.io/cran/SpadeR/man/SimilarityPair.html))
– which generates a more comprehensive set of measures and also standard
errors/confidence intervals – the main advantage of `dissCqN` is it’s
simplicity and speed, when only the original empirical C<sub>*qN*</sub>
measures are required (and also if dissimilarity is preferred to
similarity). Everything can be accomplished with a single function,
[`dissCqN()`](https://murphymv.github.io/dissCqN/reference/dissCqN.html),
which takes a matrix of assemblages x species as it’s first argument (or
a list of species interaction matrices, for network dissimilarity).

## Installation

You can install the released version of `dissCqN` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dissCqN")
```

And the development version from [GitHub](https://github.com/) with:

``` r
devtools::install_github("murphymv/dissCqN@dev")
```

## Examples

See the following vignette for a demonstration:

-   [Calculating multiple and pairwise assemblage dissimilarity for tree
    species in a tropical
    forest](https://murphymv.github.io/dissCqN/articles/dissCqN.html)

## References

Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H., & Chazdon, R. L.
(2008). A Two-Stage Probabilistic Approach to Multiple-Community
Similarity Indices. *Biometrics*, *64*(4), 1178–1186.
<https://doi.org/10/fcvn63>

Jost, L., Chao, A., & Chazdon, R. L. (2011). Compositional similarity
and beta diversity. In A. E. Magurran & B. J. McGill (Eds.), *Biological
Diversity: Frontiers in Measurement and Assessment* (pp. 66–84). Oxford
University Press.
