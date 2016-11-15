#include "basic_util.h"

#ifndef SIR_LOG_LNA_H_
#define SIR_Log_LNA_H_
using namespace Rcpp;
using namespace arma;

arma::mat Fm_log_LNA(double X,double Y,double theta1,double theta2);
List IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2);
arma::vec SIR_ODE(double X,double Y,double theta1,double theta2);

#endif /* SIR_LOG_LNA_H_*/

