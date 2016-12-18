#include "basic_util.h"
#include "SIR_LNA.h"

#define pi 3.14159265358979323846264338327950288

arma::mat SIRS_Fm_LNA_period(double X,double Y, double Z, double theta1,double theta2, double theta3,double alpha, double t){
  arma::mat M;
  double th1 = theta1 * (1 + alpha * sin(t / 40 * 2 * pi));
  M << -th1 * Y << - th1 * X << theta3 <<arma::endr
    <<th1 * Y<< th1*X - theta2<< 0 <<arma::endr
    <<0 << theta2 << - theta3 <<endr;
  return M;
}


//[[Rcpp::export()]]
List SIRS_IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2,double theta3,double alpha){
  arma::mat Sigma,F(3,3);
  Sigma.zeros(3,3);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(3);
  A << -1 << 1 << 0 <<arma::endr<<0<<-1<<1 <<arma::endr
    << 1 <<0 << -1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0.zeros(3,3);
  double t;
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    double th1 = theta1 * (1 + alpha * sin(t / 40 * 2 * pi));
    h(0)= th1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = theta2 * Traj_par(i,2);
    h(2) = theta3 * Traj_par(i,3);
    H = diagmat(h);
    F = SIRS_Fm_LNA_period(Traj_par(i,1),Traj_par(i,2), Traj_par(i,3), theta1, theta2, theta3,alpha,t);
    F0 = F0 + F*dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  A.t() * H * A ) * dt;
  }
  List Res;
  Res["expF"] = arma::expmat(F0);
  Res["Simga"] = (Sigma + 0.00000000001 * eye(3,3));
  if(Sigma.has_nan()){
    Rcout<<Traj_par<<endl;
  }
  return Res;
}


//[[Rcpp::export()]]
arma::vec SIRS_ODE(arma::vec states, arma::vec param,double t){
  double dx, dy, dz;
  double th1 = param[0] * (1 + param[3] * sin(2 * pi * t / 40.0));
  //Rcout<< th1 <<endl;
  dx = - th1 * states[0] * states[1] + param[2] * states[2];
  dy = th1 * states[0] * states[1] - param[1] * states[1];
  dz = param[1] * states[1] - param[2] * states[2];
  arma::vec res(3);
  res(0) = dx;
  res(1) = dy;
  res(2) = dz;
  //  res(2) = dz;
  return res;
}


//[[Rcpp::export()]]
List SIRS_KOM_Filter(arma::mat OdeTraj, arma::vec param,int gridsize){
  int n = OdeTraj.n_rows;
  double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  int p = OdeTraj.n_cols - 1;
  arma::cube Acube(p,p,k);
  arma::cube Scube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = SIRS_IntSigma(Traj_part,dt,param[0],param[1],param[2],param[3]);
   // Rcout<< tempres<<endl;
    Acube.slice(i) = as<arma::mat>(tempres[0]);
    Scube.slice(i) = as<arma::mat>(tempres[1]);
  }
  List Res;
  Res["A"] = Acube;
  Res["Sigma"] = Scube;
  return Res;
}

// The deterministic part of SIR after Ito transformation



//[[Rcpp::export()]]
arma::mat ODE2(arma::vec initial, arma::vec t, arma::vec param){
// XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  int n = t.n_rows;
  int p = initial.size();
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,p + 1);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,p) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  for(int i = 1; i < n; i++){
    X0 = X1;
    k1 = SIRS_ODE(X0,param,t[i-1]);
    k2 = SIRS_ODE(X0 + k1 * dt / 2, param, t[i-1]);
    k3 = SIRS_ODE(X0 + k2 * dt / 2, param, t[i-1]);
    k4 = SIRS_ODE(X0 + k3 * dt / 2, param, t[i-1]);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    OdeTraj.submat(i,1,i,p) = X1.t();
  }
  return OdeTraj;
}


//[[Rcpp::export()]]
List Traj_sim_SIRS(arma::vec initial, arma::mat OdeTraj, List Filter,double t_correct = 90){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  int k = OdeTraj.n_rows - 1;
  int p = OdeTraj.n_cols - 1;
  arma::vec X0,X1 = initial,eta0(p),eta1 = initial;

  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = initial.t();
  //  Rcout<<"test1"<<endl;
  for(int i = 0; i< k; i++){

    arma::mat Sig = Scube.slice((i));
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, p).t();
    //    Rcout<<"test2"<<endl;
    arma::mat noise = mvrnormArma(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj(i+1,0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * arma::inv(Sig) * noise;
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(det(Sig)) + l1(0,0);
    }
    //    Rcout<<"test3"<<endl;
    LNA_traj(i+1,0) = OdeTraj((i+1), 0);
    LNA_traj.submat(i+1,1,i+1,p) = X1.t();
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["loglike"] = loglike;
  return Res;
}
