# LNAphyloDyn: 
Package of the research paper: [Fitting stochastic epidemic models to gene genealogies using linear noise approximation](https://arxiv.org/abs/1902.08877)

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

### See Errors in installation?
If you see errors when installing the package on your laptop, you can simply download the package, go to the directory and run
```r
sourceCpp('src/LNA_functional.cpp')
source('R/basic.R', echo = TRUE)
source('R/coal_simulation.R', echo=TRUE)
source('R/general_change_point.R', echo=TRUE)
source('R/LNA_chpt_slice.R', echo=TRUE)
```
The above files contains most of the essential functions in the packages. 

## Vignettes
Here are some Vignettes that helps you get familiar with our work. Just download and have fun! 

1. [ConstantR0](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/constant_sim.Rmd) A short example showing how to simulate SIR trajectory and genealogy based on constant reproduction number R0, illustrating methodology ...

2. [ChangeR0](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/Changpoint_sim.Rmd) A short example showing how to simulate SIR trajectory and genealogy based time-varying reproduction number R0, illustrating methodology ...

3. [Ebola](https://github.com/MingweiWilliamTang/LNAphyloDyn/blob/master/vignettes/Ebola_sierra_leone2014.Rmd) Case studies analyzing Ebola genealogies from Sierra Leone and Liberia in the 2014 ebola outbreaks. 
