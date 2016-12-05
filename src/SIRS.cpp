#include "basic_util.h"
#include "SIR_LNA.h"

arma::mat SIRS_Fm_LNA(double X,double Y, double Z,double theta1,double theta2, double theta3){
  arma::mat M;
  M << -theta1 * Y << - theta1 * X << theta3 <<arma::endr
    <<theta1 * Y<< theta1*X - theta2<< 0 <<arma::endr
 <<0 << theta2 << - theta3 <<endr;
  return M;
}



List SIRS_IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2,double theta3){
  arma::mat Sigma,F(3,3);
  Sigma.zeros(3,3);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(3);
  A << -1 << 1 << 0 <<arma::endr<<0<<-1<<1 <<arma::endr
<< 1 <<0 << -1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0.zeros(3,3);
  for(int i = 0; i < k; i++){
    h(0)= theta1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = theta2 * Traj_par(i,2);
    h(2) = theta3 * Traj_par(i,3);
    H = diagmat(h);
    F = Fm_LNA(Traj_par(i,1),Traj_par(i,2), Traj_par(i,3), theta1, theta2, theta3);
    F0 = F0 + F*dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  A.t() * H * A ) * dt;
  }
  List Res;
  Res["expF"] = expM(F0);
  Res["Simga"] = (Sigma + 0.00000000001 * eye(3,3));
  if(Sigma.has_nan()){
    Rcout<<Traj_par<<endl;
  }
  return Res;
}


//[[Rcpp::export()]]
arma::vec SIRS_ODE(double X,double Y, double Z, double theta1,double theta2, double theta3){
  double dx, dy, dz;
  dx = - theta1*X*Y;
  dy = theta1 * X * Y - theta2 * Y;
  dz = theta2 * Y - theta3 * Z;
  arma::vec res(3);
  res(0) = dx;
  res(1) = dy;
  res(3) = dz;
  //  res(2) = dz;
  return res;
}

