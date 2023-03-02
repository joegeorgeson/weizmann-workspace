# setup on wexac
The easiest way to install and not destroy your environment is to build inside a conda
```
conda create --name brainstorm -c conda-forge r-base
conda activate brainstorm
conda install -c bioconda bioconductor-genomicranges
conda install -c bioconda bioconductor-biostrings
conda install -c bioconda bioconductor-iranges
#conda install -c bioconda bioconductor-shortread #already installed
#conda install -c bioconda bioconductor-rsamtools #already installed
#conda install -c bioconda bioconductor-s4vectors #already installed
#conda install -c bioconda bioconductor-genomicalignments #already installed

conda install -c bioconda bioconductor-plyranges #clashed with previous installs

conda install -c conda-forge r-devtools
conda install -c conda-forge r-plyr
conda install -c conda-forge r-optparse

R
devtools::install("/home/labs/schwartzlab/Collaboration/programs/brainSTORM")
# don't upgrade any package
```


# brainSTORM <img src='man/figures/logo.png' align="right" height="150" /></a>

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.0.7.1-blue.svg)](https://github.com/SchwartzLab/brainSTORM)
<!-- badges: end -->

## Description

**brainSTORM** is a package that processes STORMseq-based data, making
use of the [**txtools**](https://github.com/AngelCampos/txtools)
package.

## Installation

You can install the development version from
[GitHub](https://github.com/SchwartzLab/brainSTORM) typing in the
following commands in the R console:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("SchwartzLab/brainSTORM", build_vignettes = FALSE)
```

If you want to build the vignette, please install with the following
command (only works while connected to WEXAC).

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("SchwartzLab/brainSTORM", build_vignettes = TRUE)
```

## Further documentation

To open the “Using brainSTORM at WEXAC” vignette run the following
command.

``` r
vignette("Using_brainSTORM_at_WEXAC", package = "brainSTORM")
```

## Current limitations

## Licence

brainSTORM is developed under the [Artistic Licence
2.0](https://opensource.org/licenses/Artistic-2.0).
