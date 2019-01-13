
<!-- README.md is generated from README.Rmd. Please edit that file -->
peas
====

The purpose of this R package is to predict offspring distributions in genetic crosses. The package contains functions for defining systems of Mendelian inheritance (defined broadly), and making predictions about the distribution of offspring genotypes and phenotypes in genetic crosses. The package can handle multiple linkage groups with autosomal, sex-chromosomal or uniparental inheritance, with arbitrary numbers of loci and alleles per locus. Multiple traits can be freely defined, potentially involving pleiotropy, polygeny and epistasis.

Example
-------

Here is a simple example for how one can define a genetic system along with a trait:

``` r
library(peas)
MendelsPeas <- newGenopheno(alleleNames = list(c('Y', 'y')))
MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
MendelsPeas
#> Genetic system comprising one linkage group:
#>   Linkage group 1: autosomal, 1 locus
#>     Alleles at locus 1: Y, y
#> Phenotypes defined for the following traits:
#>   colour (trait values: yellow, green)
```

Once defined, it is easy to predict offspring genotype and phenotype distributions in crosses:

``` r
predictCross(MendelsPeas, 'Yy', 'Yy')
#> $genotypes
#>    fraction
#> YY     0.25
#> Yy     0.50
#> yy     0.25
#> 
#> $phenotypes
#>   colour fraction
#> 1  green     0.25
#> 2 yellow     0.75
```

A detailed description of the package with more examples can be found in the package vignette (after installing the package, call `browseVignettes("peas")`).

Installation
------------

To install the peas package, you can run the following two lines of code in R:

``` r
install.packages("devtools")
devtools::install_github("JanEngelstaedter/peas", build_vignettes = TRUE)
```
