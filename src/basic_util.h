#ifndef BASIC_UTIL_H_
#define BASIC_UTIL_H_

#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
using namespace Rcpp;
using namespace arma;
arma::mat inv2(arma::mat a);
arma::mat chols(arma::mat S);
arma::mat mvrnormArma(int n, arma::mat sigma);

arma::mat expM(arma::mat A);

//arma::mat ODE(arma::vec initial, arma::vec t, arma::vec param);

#endif/* BASIC_UTIL_H_ */
