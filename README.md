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
## Vignettes
1. [ConstantR0](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/constant_sim.Rmd) A short example showing how to simulate SIR trajectory and genealogy based on constant reproduction number R0, illustrating methodology ...

2. [ChangeR0](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/Changpoint_sim.Rmd) A short example showing how to simulate SIR trajectory and genealogy based time-varying reproduction number R0, illustrating methodology ...

3. [Ebola](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/Ebola_sierra_leone2014.Rmd) Case studies analyzing Ebola genealogies from Sierra Leone and Liberia in the 2014 ebola outbreaks. 
