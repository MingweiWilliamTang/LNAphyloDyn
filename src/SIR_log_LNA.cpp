#include "basic_util.h"
#include "SIR_log_LNA.h"

arma::mat Fm_log_LNA(double X,double Y,double theta1,double theta2){
  arma::mat M;
  M << theta1 * exp(Y-X)/2 << -theta1*exp(Y)-theta1/2*exp(Y-X) <<arma::endr
    <<theta1*exp(X)-theta1*exp(X-Y)/2 << theta1*exp(X-Y)/2+theta2/2*exp(-Y)<<arma::endr;
  return M;
}




List IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2){
  arma::mat Sigma,F(2,2);
  Sigma.zeros(2,2);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(2);
  A << -1 << 1 <<arma::endr<<0<<-1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0.zeros(2,2);
  for(int i = 0; i < k; i++){
    h(0)= theta1*exp(Traj_par(i,1)+Traj_par(i,2));
    h(1) = theta2*exp(Traj_par(i,2));
    H = diagmat(h);
    F = Fm_log_LNA(Traj_par(i,1),Traj_par(i,2),theta1,theta2);
    arma::vec invec(2);
    invec(0) = exp(-Traj_par(i,1));
    invec(1) = exp(-Traj_par(i,2));
    Xinv = diagmat(invec);
    F0 = F0 + F*dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma + Xinv * A.t() * H * A * Xinv) * dt;
  }
  List Res;
  Res["expF"] = expmat(F0);
  Res["Simga"] = Sigma + 0.000000001 * eye(2,2);
  return Res;
}


//[[Rcpp::export()]]
arma::vec SIR_ODE(double X,double Y,double theta1,double theta2){
  double dx,dy,dz;
  dx = -theta1 * exp(Y) - theta1*exp(Y-X)/2;
  dy = theta1 * exp(X) - theta2 - theta1 * exp(X-Y)/2 - theta2 *exp(-Y)/2;
  //  dz = theta2*exp(Y-Z) - theta2*exp(Y-2*Z)/2;
  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  //  res(2) = dz;
  return res;
}
