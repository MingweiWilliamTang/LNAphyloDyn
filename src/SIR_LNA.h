#include "basic_util.h"

#ifndef SIR_LNA_H_
#define SIR_LNA_H_
using namespace Rcpp;
using namespace arma;

arma::mat Fm_LNA(double X,double Y, double theta1,double theta2);
List IntSigma2(arma::mat Traj_par,double dt,double theta1,double theta2);
arma::mat LogTraj(arma::mat Traj);
arma::vec SIR_ODE2(double X,double Y,double theta1,double theta2);


#endif /* SIR_LNA_H_ */
