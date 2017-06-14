//@copyright Mingwei
#include "basic_util.h"

//[[Rcpp::export()]]
arma::mat inv2(arma::mat a){
  arma::mat res(2,2);
  res(0,0) = a(1,1);
  res(1,1) = a(0,0);
  res(1,0) = -a(1,0);
  res(0,1) = -a(0,1);
  double D = res(0,0) * res(1,1) - res(1,0) * res(0,1);
  return res * (1/D);
}

//[[Rcpp::export()]]
arma::mat chols(arma::mat S){
  arma::mat res(2,2);
  res(0,0) = sqrt(S(0,0));
  res(0,1) = S(1,0) / res(0,0);
  res(1,0) = 0;
  res(1,1) = sqrt(S(1,1) - res(1,0) * res(1,0));
  return res;
}

//[[Rcpp::export()]]
arma::mat mvrnormArma(int n, arma::mat sigma) {
  int ncols = sigma.n_cols;
 arma::mat Y = randn(ncols,1);
// arma::mat res = arma::chol(sigma+0.00000000001 * arma::diagmat(ones(3))) * Y;
arma::mat res;
if(n == 2){
  res = chols(sigma) * Y;
}else{
  res = arma::chol(sigma + 0.000000001 * arma::diagmat(ones(3)) ) * Y;
}
  if(res.has_nan()){
    Rcout<<"666"<<sigma<<endl;
  }
  return res;
}

//[[Rcpp::export()]]
arma::mat mvrnormArma2(int n, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = randn(ncols,1);
  arma::mat res = arma::chol(sigma + 0.000000001 * arma::diagmat(ones(3)) ) * Y;
//arma::mat res = chols(sigma) * Y;
  if(res.has_nan()){
    Rcout<<"666"<<sigma<<endl;
  }
  return res;
}



//[[Rcpp::export()]]
arma::mat expM(arma::mat A){
  arma::mat M;
  double c = exp((A(0,0) + A(1,1)) / 2);
  double d = (A(0,0) - A(1,1)) * (A(0,0) - A(1,1)) + 4 * A(1,0) * A(0,1);
  double delta = sqrt(fabs(d));
  if(d >= 0.000000000001){
    double m11 = c * ((delta * cosh(delta / 2)) + (A(0,0) - A(1,1)) * sinh(delta / 2));
    double m12 = 2 * A(0,1) * c * sinh(delta / 2);
    double m21 = 2 * A(1,0) * c * sinh(delta / 2);
    double m22 = c * ((delta * cosh(delta / 2)) + (A(1,1) - A(0,0)) * sinh(delta / 2));
    M<< m11/delta << m12/delta <<arma::endr << m21/delta << m22/delta <<arma::endr;
  }else if(d < -0.0000000001){
    double m11 = c * ((delta * cos(delta / 2)) + (A(0,0) - A(1,1)) * sin(delta / 2));
    double m12 = 2 * A(0,1) * c * sin(delta / 2);
    double m21 = 2 * A(1,0) * c * sin(delta / 2);
    double m22 = c * ((delta * cos(delta / 2)) + (A(1,1) - A(0,0)) * sin(delta / 2));
    M<< m11/delta << m12/delta <<arma::endr << m21/delta << m22/delta <<arma::endr;
  }else{
    double m11 = c * (1 + (A(0,0) - A(1,1)) /2 );
    double m12 = c * A(0,1);
    double m21 = c * A(1,0);
    double m22 = c * (1 - (A(0,0) - A(1,1)) /2 );
    M<< m11 << m12 <<arma::endr << m21 << m22 <<arma::endr;
  }
  return M;
}

