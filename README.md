# miaViz

<!-- badges: start -->

[![R-CMD-check-Bioc-devel](https://github.com/microbiome/miaViz/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/microbiome/miaViz/actions)
[![Codecov test
coverage](https://codecov.io/gh/microbiome/miaViz/branch/master/graph/badge.svg)](https://codecov.io/gh/microbiome/miaViz?branch=master)

<!-- badges: end -->

## Microbiome Analysis Plotting and Visualization

The scope of this package is the plotting and visualization of microbiome data.
The main classes for interfacing is the `TreeSummarizedExperiment` class.

The [package homepage](https://microbiome.github.io/miaViz/) provides
some online tutorials and examples.

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

# Code of conduct

Please note that the miaViz project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
