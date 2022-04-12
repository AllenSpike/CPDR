
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cancer Personalized Drug Recommendation (CPDR)

<!-- badges: start -->
<!-- badges: end -->

A recommendation tool of personalized drugs for patients with cancer by
reversing individual disease signatures

## Installation

Before library installation install required Bioconductor and CRAN
packages through this code:

``` r
bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5')

#For R version 3.5> use BiocManager to install required bioconductor packages: 
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

packages=c('magrittr','dplyr','ggplot2','doParallel','foreach','Rfast','data.table')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
```

Now you can install the development version of CPDR with:

``` r
install.packages('devtools')
devtools::install_github("AllenSpike/CPDR")
```

## Additional data

By default, CPDR uses LINCS\_978 from octad.db. To acquire reliable
results, It is suggested to download the LINCS\_all\_processed data
(6.5Gb) manually from <https://zenodo.org/record/5880026#.YgZos9--uUl>.

## Example

The vignette shows the usecase of screening effective compounds
targeting individual patient with colorectal cancer.
