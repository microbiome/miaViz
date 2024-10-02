# miaViz <img src="man/figures/mia_logo.png" align="right" width="120" />

<!-- badges: start -->
[![Platforms](http://bioconductor.org/shields/availability/release/miaViz.svg)](https://bioconductor.org/packages/release/bioc/html/miaViz.html)
[![rworkflows](https://github.com/microbiome/miaViz/actions/workflows/rworkflows.yml/badge.svg?branch=devel)](https://github.com/microbiome/miaViz/actions)
[![Bioc-release](http://bioconductor.org/shields/build/release/bioc/miaViz.svg)](http://bioconductor.org/packages/release/bioc/html/miaViz.html)
[![Bioc-age](http://bioconductor.org/shields/years-in-bioc/miaViz.svg)](https://bioconductor.org/packages/release/bioc/html/miaViz.html#since)
[![Codecov test
coverage](https://codecov.io/gh/microbiome/miaViz/branch/devel/graph/badge.svg)](https://codecov.io/gh/microbiome/miaViz?branch=devel)
[![Dependencies](http://bioconductor.org//shields/dependencies/release/miaViz.svg)](https://bioconductor.org/packages/release/bioc/html/miaViz.html)
<!-- badges: end -->

## Microbiome Analysis Plotting and Visualization

The scope of this package is the plotting and visualization of microbiome data.
The main classes for interfacing is the `TreeSummarizedExperiment` class.

## Using the package

Online tutorials and examples are available at:

- [Package homepage](https://microbiome.github.io/miaViz/) 
- [Orchestrating microbiome analysis with R/Bioconductor online book](https://microbiome.github.io/OMA)


## Contribution

Feel free to contribute by forking and opening a pull request. Please make sure
that required data wrangling should be designed as reusable as possible and
potentially find a better home in the [`mia`](https://github.com/FelixErnst/mia)
package.

Additionally, please make sure that working examples are included and that 
vignetted make use of added functions in either `miaViz` or the
[`TreeSummarizedExperiment`](https://github.com/fionarhuang/TreeSummarizedExperiment)
package.

## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `release` for release. 
(https://guides.github.com/introduction/flow/)

## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miaViz")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("miaViz")
```

# Code of conduct

Please note that the miaViz project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
