Overview
--------

SIGN fasciliotates classification and clustering of biological samples relyign on expression pattersn of biological pathways. A new measure of pathway expression pattern similarity (TSC) was introduced in the package. The package has been developed and tested for RNA-seq profiles of the cells. However, it can be used for other sequencig profiles with continuous values for each feature (gene, protein, cis-regulatory elements, etc.)


Installation
------------

``` r
# Install from CRAN
install.packages('SIGN')

# Installing the development version from GitHub:
# install.packages("devtools")
devtools::install_github("bhklab/SIGN")
