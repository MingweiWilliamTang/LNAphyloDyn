#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

using namespace Rcpp;
using namespace arma;

typedef arma::vec (*parat)(double, arma::vec, arma::vec, arma::ivec);
typedef arma::vec (*ODE_one)(arma::vec,double,arma::vec,arma::vec, arma::ivec, std::string);
typedef arma::mat (*F_fun)(arma::vec, arma::vec, std::string);
typedef arma::vec (*h_fun)(arma::vec, arma::vec, std::string);
//typedef List (*SigmaInt)(arma::mat,arma::vec,arma::vec, arma::ivec, std::string, std::string, std::string);

//[[Rcpp::export()]]
arma::vec param_transform(double t, arma::vec param, arma::vec x_r, arma::ivec x_i){
  /**
  * x_r = (N,cht1, cht2, ...)
  * theta = (R0,gamma, lambda, ch1,ch2,...)
  * x_i = (nch,nparam)
  */
  arma::vec param2 = param;
  double R0 = param[0];
  int i = 0;
  int nch = x_i[0];
  if(nch > 0){
    while(x_r[i + 1] <= t){
      R0 *= param[x_i[1] + i];
      if(i == nch - 1) break;
      i ++;
    }
    param2[0] = R0 * param[1] / x_r[0];
  }
  return param2;
}

//[[Rcpp::export()]]
XPtr<parat> transformPtr(std::string trans = "changepoint"){
    if(trans == "changepoint"){
      return(XPtr<parat>(new parat(&param_transform)));
    }else{
      return XPtr<parat>(R_NilValue); // runtime error as NULL no XPtr
    }
}




//[[Rcpp::export()]]
arma::vec ODE_SIR_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i,std::string transP = "changepoint", std::string transX = "standard"){
  double dx, dy;
  arma::vec res(2);
  XPtr<parat> param_trans = transformPtr(transP);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    dx = - thetas[0] * states[0] * states[1];
    dy = - thetas[0] * states[0] * states[1] - thetas[1] * states[1];
  }else if(transX == "log"){
    dx = -thetas[0] * exp(Y) - thetas[0]*exp(Y-X)/2;
    dy = thetas[0] * exp(X) - thetas[1] - thetas[0] * exp(X-Y)/2 - thetas[1] *exp(-Y)/2;
  }
  res(0) = dx;
  res(1) = dy;
  return res;
}




//[[Rcpp::export()]]
arma::vec ODE_SIRS_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i){
  double dx, dy;
  XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  dx = - thetas[0] * states[0] * states[1];
  dy = - thetas[0] * states[0] * states[1] - thetas[1] * states[1];

  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}


//[[Rcpp::export()]]
arma::mat SIR_F(arma::vec states,arma::vec thetas,std::string transX){
  arma::mat M;
  arma::mat A;
  arma::vec h;
  A << -1 << 1  <<arma::endr<<0<< -1 << arma::endr;
  //XPtr<parat> param_trans = transformPtr();
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  //arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
  if(transX == "standard"){
    M << - thetas[0] * states[1] << - thetas[0] * states[0] << arma::endr
      << thetas[0] * states[1] << thetas[0] * states[0] - thetas[1] << arma::endr;
  }else if(transX == "log"){
    M << thetas[0] * exp(states(1)-states(0))/2 << -thetas[0]*exp(states(1))-thetas[0]/2*exp(states(1)-states(0)) <<arma::endr
      <<thetas[0]*exp(states(0))-thetas[0]*exp(states(0)-states(1))/2 << thetas[0]*exp(states(0)-states(1))/2+thetas[1]/2*exp(-states(1))<<arma::endr;
  }
  return M;
}


//[[Rcpp::export()]]
arma::vec SIR_h(arma::vec states,arma::vec thetas,std::string transX = "standard"){
  arma::vec h(2);
  if(transX == "standard"){
    //XPtr<parat> param_trans = transformPtr();
    //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
    //arma::vec thetas = (*param_trans)(tstates(0),param,x_r,x_i);
    h(0)= thetas[0] * states(0) * states(1);
    h(1) = thetas[1] * states(2);
  }else if(transX == "log"){
    h(0)= thetas[0] * exp(states(1)) * exp(states(2));
    h(1) = thetas[1] * exp(states(2));
  }
  return h;

  }



//[[Rcpp::export()]]
arma::vec ODE_SIR_log_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
                           std::string transP = "changepoint"){
  double dx, dy;
  XPtr<parat> param_trans = transformPtr(transP);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  arma::vec thetas = (*param_trans)(t,param,x_r,x_i);

  dx = - thetas[0] * states[0] * states[1];
  dy = - thetas[0] * states[0] * states[1] - thetas[1] * states[1];

  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}
/*
arma::mat F_SIR_log(arma::vec states, arma::vec transParam){
  arma::mat M;
  double X = states[0],Y = states[1];
  M << transParam[0] * exp(Y-X)/2 << -transParam[0]*exp(Y)-transParam[0]/2*exp(Y-X) <<arma::endr
    <<transParam[0]*exp(X)-transParam[0]*exp(X-Y)/2 << transParam[0]*exp(X-Y)/2+transParam[1]/2*exp(-Y)<<arma::endr;
  return M;
}
*/



//[[Rcpp::export()]]
XPtr<F_fun> F_funPtr(std::string model = "SIR"){
  if(model == "SIR" ){
    return XPtr<F_fun>(new F_fun(&SIR_F));
  }else if(model == "SIRS"){
    return XPtr<F_fun>(new F_fun(&SIR_F));
  }else{
    return XPtr<F_fun>(R_NilValue); // runtime error as NULL no XPtr
  }
  }

//[[Rcpp::export()]]
XPtr<h_fun> h_fun_Ptr(std::string model = "SIR"){
  if(model == "SIR"){
    return(XPtr<h_fun>(new h_fun(&SIR_h)));
  }else{
    return XPtr<h_fun>(R_NilValue); // runtime error as NULL no XPtr
  }
}


//[[Rcpp::export()]]
arma::mat ODE_rk45(arma::vec initial, arma::vec t, arma::vec param,
                           arma::vec x_r, arma::ivec x_i, std::string model){
  // XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  //double N = initial[0] + initial[1] + initial[2];
  int n = t.n_rows;
  int p = initial.size();
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,p + 1);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,p) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  for(int i = 1; i < n; i++){
    X0 = X1;
    k1 = ODE_one2(X0,param,t[i-1], x_r,x_i,&param_transform);
    k2 = ODE_one2(X0 + k1 * dt / 2,param, t[i-1], x_r,x_i,param_transform);
    k3 = ODE_one2(X0 + k2 * dt / 2,param, t[i-1], x_r,x_i,param_transform);
    k4 = ODE_one2(X0 + k3 * dt / 2,param, t[i-1], x_r,x_i,param_transform);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    OdeTraj.submat(i,1,i,p) = X1.t();
  }
  return OdeTraj;
}


//[[Rcpp::export()]]
List SigmaF(arma::mat Traj_par,arma::vec param,
                  arma::vec x_r, arma::ivec x_i,std::string model,std::string transX, std::string transP){
  int p = Traj_par.n_cols - 1;
  arma::mat Q = arma::eye(p,p);

  arma::mat Sigma,F,F0 = arma::eye(p,p);
  Sigma.zeros(p,p);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h,state;
  if(model == "SIR"){
    A << -1 << 1  <<arma::endr<<0<< -1 << arma::endr;
  }else if(model == "SEIR"){
    A << -1 << 0 << 1 <<arma::endr << 0 << 1 << -1 <<arma::endr << 0 << -1 << 0 << endr;
  }
  arma::mat H;
  //F0 << 1 <<0 <<arma::endr << 0 << 1 <<endr;
  double t;
  double dt = Traj_par(1,0) - Traj_par(0,0);
  XPtr<parat> param_trans = transformPtr(transP);
  XPtr<h_fun> h_Function = h_fun_Ptr(model);
  XPtr<F_fun> F_function = F_funPtr(model);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    arma::vec thetas = (*param_trans)(t,param,x_r,x_i);
    state = Traj_par.row(i).subvec(1,p);
    if(transX == "log"){
      for(int i = 0; i < p; i ++){
        Q(i,i) = exp(-state(i));
      }
    }
    h = (*h_Function)(state, thetas, transX);
    H = diagmat(h);
    F = (*F_function)(state, thetas, transX);
    F0 = F0 + (F * F0) * dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  Q * A.t() * H * A * Q ) * dt;
  }
  List Res;
  Res["expF"] = F0;
  Res["Simga"] = Sigma;
  if(Sigma.has_nan()){
    Rcout<<param<<endl;
    Rcout<<Traj_par.row(0)<<endl;
  }
  return Res;
}
