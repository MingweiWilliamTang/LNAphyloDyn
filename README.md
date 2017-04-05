# LNAphyloDyn: 
Package of the research: Bayesian inference on coalescent-based phylodynamics

## Installation
### Prerequisites
You need to install these R packages before you install LNAPhlyodyn package: abind, smfsb, Rcpp, RcppArmadillo. You can easily install these packages by running the following code in your R console. 

```r
install.packages("abind")
install.packages("smfsb")
install.packages("Rcpp")
install.packages("RcppArmadillo")
```

### Install LNAPhylodyn package
After installing the prerequisites packages, LNAphyloDyn package can be installed by running the following code in your console

```r
Sys.setenv("PKG_LIBS" = "-llapack")
devtools::install_github("MingweiWilliamTang/LNAPhyloDyn")
```
