#include "basic_util.h"
#include "SIR_LNA.h"

arma::mat Fm_LNA(double X,double Y, double theta1,double theta2){
  arma::mat M;
  M << -theta1 * Y << - theta1 * X <<arma::endr
    <<theta1 * Y<< theta1*X - theta2<<arma::endr;
  return M;
}



List IntSigma2(arma::mat Traj_par,double dt,double theta1,double theta2){
  arma::mat Sigma,F(2,2);
  Sigma.zeros(2,2);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(2);
  A << -1 << 1 <<arma::endr<<0<<-1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0.zeros(2,2);
  for(int i = 0; i < k; i++){
    h(0)= theta1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = theta2 * Traj_par(i,2);
    H = diagmat(h);
    F = Fm_LNA(Traj_par(i,1),Traj_par(i,2),theta1,theta2);
    F0 = F0 + F*dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  A.t() * H * A ) * dt;
  }
  List Res;
  Res["expF"] = expM(F0);
  Res["Simga"] = (Sigma + 0.00000000001 * eye(2,2));
  if(Sigma.has_nan()){
    Rcout<<Traj_par<<endl;
  }
  return Res;
}

//[[Rcpp::export()]]
arma::mat LogTraj(arma::mat Traj){
  int p = Traj.n_cols - 1;
  arma::mat logTraj(Traj.n_rows, Traj.n_cols);
  logTraj.col(0) = Traj.col(0);
  logTraj.cols(1,p) = log(Traj.cols(1,p));
  return logTraj;
}

//[[Rcpp::export()]]
arma::vec SIR_ODE2(double X,double Y,double theta1,double theta2){
  double dx,dy;
  dx = - theta1*X*Y;
  dy = theta1 * X * Y - theta2 * Y;
  //  dz = theta2*exp(Y-Z) - theta2*exp(Y-2*Z)/2;
  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  //  res(2) = dz;
  return res;
}

