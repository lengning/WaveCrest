# WaveCrest


WaveCrest is a statistical approach to
reconstruct gene expression trajectory in single cell RNA-seq experiments with ordered conditions.
WaveCrest contains two modules - the first module implements an extended nearest insertion (ENI) algorithm that
searches for optimal cell orders, and the second module implements a spline fitting module
that can be used to identify additional dynamic genes.

## Installation

SCPattern R package and its vignette could be found at https://github.com/lengning/SCPattern/tree/master/package

To install SCPattern, in R run: 

> install.packages("devtools")

> library(devtools)

> install_github("lengning/SCPattern/package/SCPattern")

Or install locally.
