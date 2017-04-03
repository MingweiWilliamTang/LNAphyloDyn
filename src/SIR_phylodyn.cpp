//@copyright Mingwei
//#include<RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
#include "basic_util.h"
#include "SIR_LNA.h"
#include "SIR_log_LNA.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
#include<iostream>
using namespace Rcpp;
using namespace arma;
# define pi           3.1415926535897932384

typedef arma::vec (*ODEfuncPtr)(double,double,double,double);
typedef List (*IntSfun)(arma::mat,double,double,double);




XPtr<ODEfuncPtr> putFunPtrInXPtr(std::string funname){
  if(funname == "log"){
    return(XPtr<ODEfuncPtr>(new ODEfuncPtr(&SIR_ODE)));
  }else{
    return(XPtr<ODEfuncPtr>(new ODEfuncPtr(&SIR_ODE2)));
  }
}


XPtr<IntSfun> putIntSFun(std::string funname){
  if(funname == "log"){
    return XPtr<IntSfun>(new IntSfun(&IntSigma));
  }else{
    return XPtr<IntSfun>(new IntSfun(&IntSigma2));
  }
}

//[[Rcpp::export()]]
arma::vec betaDyn(double beta, double alpha, arma::vec times, double period = 40){
   arma::vec betat(times.n_elem);
  for(int i = 0; i < times.n_elem; i ++){
   betat[i] = beta * (1 + alpha * sin(2 * pi * times[i] / period));
  }
  return betat;
}


//[[Rcpp::export()]]
arma::mat SIRS2_period_SDE(arma::vec init, double N, arma::vec param, arma::vec t, double period = 40){
  int n = t.size();
  int p = init.size();
  arma::mat H;
  arma::mat Traj(n,p + 1);
  Traj.col(0) = t;
  Traj.submat(0,1,0,p) = init.t();
  arma::mat A;
  A << -1 << 0 << 1 << 0 << endr
    << 0 << -1 << 1 << 0 << endr
    << -1 << 1 << 0 << 0 << endr
    << 0 << 0 << -1 << 1 << endr
    << 1 << 0 << 0 << -1 << endr;
  double dt  = t(1) - t(0);
  double beta1 = param(0);
  double beta2 = param(1);
  double beta3 = param(2);
  double gamma = param(3);
  double mu = param(4);
  double alpha = param(5);
  for( int i = 0; i < n-1; i ++){
    // double S_new, I_new, R_new;
    double betat1 = beta1 * (1 + alpha * (sin(2 * pi * t(i) / period)));
    double betat2 = beta2 * (1 + alpha * (sin(2 * pi * t(i) / period)));
    arma::mat noises = randn(5,1) * sqrt(dt);
    arma::vec h(5);
    h(0) = betat1 * Traj(i,1) * Traj(i,3);
    h(1) = betat2 * Traj(i,2) * Traj(i,3);
    h(2) = beta3 * Traj(i,1);
    h(3) = gamma * Traj(i,3);
    h(4) = mu * Traj(i,4);
    //Rcout<<h<<endl;
    arma::mat mean = A.t() * h * dt;
    //mean = mean + A.t() * arma::diagmat(arma::sqrt(h)) * noises;
    Traj.submat(i+1,1,i+1,p) = Traj.submat(i,1,i,p) + mean.t();
  }
  return Traj;
}



//[[Rcpp::export()]]
double log_like_traj(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                     int gridsize,double t_correct = 90){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows - 1;
  int p = SdeTraj.n_cols - 1;
  double loglik = 0;
//  int id1,id0;
  arma:vec Xd1, Xd0;

   Xd0 = Xd1 = (SdeTraj.submat(0,1,0,p)-OdeTraj.submat(0,1,0,p)).t();
// Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
  //  id0 = i * gridsize;
  //  id1 = (i+1) * gridsize - 1;

    arma::mat SigInv = arma::inv(Scube.slice(i) );
    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),p)-OdeTraj.submat(i+1,1,i+1,p)).t();

    if(SdeTraj(i+1,0) <= t_correct){
    arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
     // Rcout <<  ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0)) <<endl;
    loglik += -log(arma::det(Scube.slice(i))) - 0.5 * INexp(0,0);
    }
/*
  arma::mat A = Acube.slice(i);
  arma::mat SigInv = inv2(Scube.slice(i)+0.00001 * eye(3,3));
  arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
  loglik += log(arma::det(SigInv)) * 1.5 - 0.5 * INexp(0,0);
 */
  }

  return loglik;
}



/*
arma::mat Fm_log_LNA(double X,double Y,double Z,double theta1,double theta2){
  arma::mat M;
  M << theta1 * exp(Y-X)/2 << -theta1*exp(Y)-theta1/2*exp(Y-X) << 0 <<arma::endr
    <<theta1*exp(X)-theta1*exp(X-Y)/2 << theta1*exp(X-Y)/2+theta2/2*exp(-Y)<< 0 <<arma::endr
    <<0<<theta2*exp(Y-Z)-theta2*exp(Y-2*Z)/2 << theta2*exp(Y-2*Z)-theta2*exp(Y-Z)<<arma::endr;
  return M;
  }
*/


/*
List IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2){
  arma::mat Sigma,F(3,3);
  Sigma.zeros(3,3);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(2);
  A << -1 << 1 <<0 <<arma::endr<<0<<-1<<1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0.zeros(3,3);
  for(int i = 0; i < k; i++){
    h(0)= theta1*exp(Traj_par(i,1)+Traj_par(i,2));
    h(1) = theta2*exp(Traj_par(i,2));
    H = diagmat(h);
    F = Fm_log_LNA(Traj_par(i,1),Traj_par(i,2),Traj_par(i,3),theta1,theta2);
    arma::vec invec(3);
    invec(0) = exp(-Traj_par(i,1));
    invec(1) = exp(-Traj_par(i,2));
    invec(2) = exp(-Traj_par(i,3));
    Xinv = diagmat(invec);
    F0 = F0 + F*dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma + Xinv * A.t() * H * A * Xinv) * dt;
  }
  List Res;
  Res["expF"] = expmat(F0);
  Res["Simga"] = Sigma;
  return Res;
}
*/

//[[Rcpp::export()]]
List SIR_log_KOM_Filter2(arma::mat OdeTraj,double theta1,double theta2,int gridsize,std::string funname = "standard"){
  XPtr<IntSfun> IntS = putIntSFun(funname);
  int n = OdeTraj.n_rows;
  double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  arma::cube Acube(2,2,k);
  arma::cube Scube(2,2,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
//    Rcout<<i<<endl;
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,2);
    List tempres((*IntS)(Traj_part,dt,theta1,theta2));
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
arma::mat ODE(arma::vec initial, arma::vec t, arma::vec param, std::string funname){
  XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  int n = t.n_rows;
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,3);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,2) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  for(int i = 1; i < n; i++){
    X0 = X1;
    k1 = (*SIR_ODEfun)(X0[0],X0[1],param[0],param[1]);
    k2 = (*SIR_ODEfun)(X0[0]+k1[0]*dt/2,X0[1]+k1[1]*dt/2,param[0],param[1]);
    k3 = (*SIR_ODEfun)(X0[0]+k2[0]*dt/2,X0[1]+k2[1]*dt/2,param[0],param[1]);
    k4 = (*SIR_ODEfun)(X0[0]+k3[0]*dt/2,X0[1]+k3[1]*dt/2,param[0],param[1]);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    OdeTraj.submat(i,1,i,2) = X1.t();
  }
  return OdeTraj;
}

//[[Rcpp::export()]]
List Traj_sim(arma::vec initial, arma::mat OdeTraj, List Filter,double t_correct = 90){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  int k = OdeTraj.n_rows-1;

  arma::vec X0,X1 = initial,eta0(2),eta1 = initial;

  double loglike = 0;
  arma::mat LNA_traj(k+1,3);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,2) = initial.t();
//  Rcout<<"test1"<<endl;
  for(int i = 0; i< k; i++){

    arma::mat Sig = Scube.slice((i));
 //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, 2).t();
  //  Rcout<<"test2"<<endl;
    arma::mat noise = mvrnormArma(2,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
   if(OdeTraj(i+1,0) <= t_correct){
    arma::mat l1 = (-0.5) * noise.t() * inv2(Sig) * noise;
 //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
  //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
    loglike += -log(det(Sig)) + l1(0,0);
   }
//    Rcout<<"test3"<<endl;
    LNA_traj(i+1,0) = OdeTraj((i+1), 0);
    LNA_traj.submat(i+1,1,i+1,2) = X1.t();
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["loglike"] = loglike;
  return Res;
}

// Simulate a LNA trajectory and give the corresponding loglikelihood at the same time

//[[Rcpp::export()]]
List Traj_sim_ez(arma::vec initial, arma::vec times,double theta1, double theta2,
        int gridsize,double t_correct = 90,std::string funname = "standard"){
  int k = times.n_rows / gridsize;
  arma::vec param(2);
  param(0) = theta1;
  param(1) = theta2;
  arma::mat OdeTraj_thin = ODE(initial,times,param,funname);
  arma::mat OdeTraj(k+1,3);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,2) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,2);
  }
  List Filter = SIR_log_KOM_Filter2(OdeTraj_thin,theta1,theta2,gridsize,funname);

  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
   k = OdeTraj.n_rows-1;

  arma::vec X0,X1 = initial,eta0(2),eta1 = initial;

  double loglike = 0;
  arma::mat LNA_traj(k+1,3);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,2) = initial.t();
  //  Rcout<<"test1"<<endl;
  for(int i = 0; i< k; i++){

    arma::mat Sig = Scube.slice((i));
    if(Sig(0,0)<0){
      Rcout<<i<<endl;
    }
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, 2).t();
    //   Rcout<<"test2"<<endl;
    arma::mat noise = mvrnormArma(2,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
   if(OdeTraj((i+1), 0) <= t_correct){
    arma::mat l1 = (-0.5) * noise.t() * inv2(Sig) * noise;
    //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
    //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
    loglike += -log(det(Sig)) + l1(0,0);
   }
    //   Rcout<<"test3"<<endl;
    LNA_traj(i+1,0) = OdeTraj((i+1), 0);
    LNA_traj.submat(i+1,1,i+1,2) = X1.t();
  }

  List Res;
  Res["SimuTraj"] = LNA_traj;
  Res["loglike"] = loglike;
  return Res;
}

/*
 * In put a SDE path trajectory of log-transformed SIR and the parameter value,
 return the loglikelihood
 */
//[[Rcpp::export()]]
double log_like_traj2(arma::mat SdeTraj,arma::vec times,arma::vec state,
                      double theta1,double theta2,int gridsize,double t_correct = 90,std::string funname = "standard"){

// generate the ODE path, calculate mean and covariance in the transition probability
  int k = SdeTraj.n_rows-1;
  arma::vec param(2);
  param(0) = theta1;
  param(1) = theta2;

  arma::mat OdeTraj_thin = ODE(state,times,param,funname);
  arma::mat OdeTraj(k+1,3);
  for(int i = 0; i< SdeTraj.n_rows; i++){
    OdeTraj.submat(i,0,i,2) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,2);
  }
  List Filter = SIR_log_KOM_Filter2(OdeTraj_thin,theta1,theta2,gridsize,funname);
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);


  double loglik = 0;
  //  int id1,id0;

  arma:vec Xd1, Xd0;
  Xd0 = Xd1 = (SdeTraj.submat(0,1,0,2)-OdeTraj.submat(0,1,0,2)).t();
  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();

  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
    //  id0 = i * gridsize;
    //  id1 = (i+1) * gridsize - 1;

    arma::mat SigInv = inv2(Scube.slice(i) );
    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),2)-OdeTraj.submat(i+1,1,i+1,2)).t();
    if(OdeTraj(i+1,0) <= t_correct){
    arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
    loglik += log(arma::det(SigInv)) - 0.5 * INexp(0,0);
    }
    /*
     arma::mat A = Acube.slice(i);
     arma::mat SigInv = inv2(Scube.slice(i)+0.00001 * eye(3,3));
     arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
     loglik += log(arma::det(SigInv)) * 1.5 - 0.5 * INexp(0,0);
     */
  }

  return loglik;
}

/*
 * Input the initial list, return the log coalescent likelihood
 *
 */
//[[Rcpp::export()]]
double coal_loglik2(List init, arma::mat f1, double t_correct, double lambda, int gridsize = 1){
  int n0 = 0;
  arma::vec ff((f1.n_rows-1) * gridsize + 1), tf((f1.n_rows-1) * gridsize + 1);
  for(int j = 0; j < (f1.n_rows - 1); j++){
    ff(j*gridsize) = f1(j,2);
    tf(j*gridsize) = f1(j,0);
    double df = (exp(f1(j+1,2)) - exp(f1(j,2)))/gridsize;
    double dt = (f1(j+1,0) - f1(j,0))/gridsize;
    for(int k = 1; k < gridsize; k ++){
      ff(j * gridsize + k) = log(exp(ff(j*gridsize)) + k * df);
      tf(j * gridsize + k) = tf(j*gridsize) + k * dt;
    }
  }
  ff((f1.n_rows-1) * gridsize) = f1(f1.n_rows-1,2);
  while(tf(n0) <= t_correct){
    n0 ++;
  }
//  Rcout<<n0<<endl;
  arma::vec f2(n0);
  for(int i = n0; i>0; i--){
    f2(n0-i) = ff(i);
  }
//  Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  if(as<int>(init[9]) != f2.n_rows){
    Rcout<<"Incorrect length for f"<<endl;
    Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);

  int start = 0;
  for(int i = 0; i< f2.n_rows; i++){
    for(int j = 0; j<gridrep(i);j++){
      f(start ++) = f2(i);
    }
  }
  arma::vec ll = -lambda * (as<vec>(init[2]) % as<vec>(init[3]) % arma::exp(-f)) +\
    as<vec>(init[4]) % (log(lambda) -f);
  return sum(ll);
}


//[[Rcpp::export()]]
double coal_loglik(List init, arma::mat f1, double t_correct, double lambda, int gridsize = 1){

  int n0 = 0;
  while(f1(n0,0) < t_correct){
    n0 ++;
  }
  //  Rcout<<n0<<endl;
  arma::vec f2(n0 + 1);
  for(int i = n0; i>0; i--){
    f2(n0-i) = f1(i,2);
  }
  //  Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  if(as<int>(init[9]) != f2.n_rows){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);

  int start = 0;
  for(int i = 0; i< f2.n_rows; i++){
    for(int j = 0; j<gridrep(i);j++){
      f(start ++) = f2(i);
    }
  }
  arma::vec ll = -lambda * (as<vec>(init[2]) % as<vec>(init[3]) % arma::exp(-f)) +\
    as<vec>(init[4]) % (log(lambda) -f);
  return sum(ll);
}



//[[Rcpp::export()]]
double volz_loglik(List init, arma::mat f1, double t_correct, double betaN, int gridsize = 1){

  int n0 = 0;
  while(f1(n0,0) < t_correct){
    n0 ++;
  }

  if(f1.submat(0,1,n0,2).min() < 0){
    return -10000000;
  }
  //  Rcout<<n0<<endl;
  // Rcout<<f1 <<endl;
  arma::mat f2(n0 + 1,2);
  for(int i = n0; i >= 0; i--){
    f2.submat((n0-i),0,(n0-i),1) = f1.submat(i,1,i,2);
  }
   // Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  if(as<int>(init[9]) != f2.n_rows){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);
  arma::vec s(k);
  int start = 0;
  for(int i = 0; i< f2.n_rows; i++){
    for(int j = 0; j<gridrep(i);j++){
      f(start) = f2(i,1);
      s(start) = f2(i,0);
      start ++;
    }
  }

  arma::vec ll = -2 * betaN * (as<vec>(init[2]) % as<vec>(init[3]) % (arma::exp(-f)) % (arma::exp(s))) +\
    as<vec>(init[4]) % (log(betaN) -f+s);
  Rcout<< (arma::exp(-f)) % (arma::exp(s))<<endl;
  return sum(ll);
}


//[[Rcpp::export()]]
double volz_loglik_nh(List init, arma::mat f1, arma::vec betaN, double t_correct, int gridsize = 1){

  int n0 = as<int>(init[9]) - 1;
 int L = 0;
  while(f1(L,0) < t_correct){
    L ++;
  }
  if(f1.submat(0,1,n0,2).min() < 0){
    return -10000000;
  }
  // Rcout<<f1 <<endl;
  arma::mat f2(n0,2);
  arma::vec betanh(n0);
  for(int i = 1; i<= n0; i++){
    f2.submat((n0-i),0,(n0-i),1) = f1.submat(L-i+1,1,L-i+1,2);
    betanh(n0-i) = betaN(L-i+1);
  }
  // Rcout<<f2.n_rows<<"\t"<<as<int>(init[9])<<endl;
  if(as<int>(init[9]) != (f2.n_rows+1)){
    Rcout<<"Incorrect length for f"<<endl;
  }

  arma::vec gridrep;
  gridrep = as<vec>(init[6]);
  int k = sum(as<arma::vec>(init[6]));

  arma::vec f(k);
  arma::vec s(k);
  arma::vec betas(k);
  int start = 0;
  for(int i = 0; i < f2.n_rows; i++){
    for(int j = 0; j < gridrep(i);j++){
      f(start) = f2(f2.n_rows-i-1,1);
      s(start) = f2(f2.n_rows-i-1,0);
      betas(start) = betanh(f2.n_rows-i-1);
      start ++;
    }
  }

  arma::vec ll = -2 * (as<vec>(init[2]) % as<vec>(init[3]) % betas % (arma::exp(-f)) % (arma::exp(s))) +\
    as<vec>(init[4]) % (log(betas) -f+s);
  //Rcout<< ll <<endl;
  return sum(ll);
}



/*
 * Input the current trajectory and the ODE trajectory as well as the
 *
 */
/*
arma::mat ESlice(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                 List init, double t_correct, double lambda=10, int reps=1, int gridsize = 100){
  // OdeTraj is the one with low resolution
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,2) - OdeTraj.cols(1,2);
    //f_cur_centered(0,1)=0;
    //f_cur_centered(0,0)=0;
    //simulate a new trajectory
    List v = Traj_sim(state,OdeTraj,FTs);
    arma::mat v_traj = as<mat>(v[0]).cols(1,2) -  OdeTraj.cols(1,2);
    if(v_traj.has_nan()){
      Rcout<<v_traj<<endl;
      Rcout<<OdeTraj<<endl;
    }
    double u = R::runif(0,1);
    logy = coal_loglik(init,f_cur,t_correct,lambda,gridsize) + log(u);

    double theta = R::runif(0,2 * pi);
    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
    newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    int i = 0;
    while(coal_loglik(init,newTraj,t_correct,lambda,gridsize) <= logy){
      // shrink the bracket
      i = 1;
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
      // newTraj.col(0) = f_cur.col(0);
      newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    }
    f_cur = newTraj;
    //Rcout<< logy - coal_loglik(init,newTraj,t_correct,lambda,gridsize)<<"1"<<endl;
  }
  return newTraj;
}

*/


//[[Rcpp::export()]]
arma::mat ESlice(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                  List init, double t_correct, double lambda=10, int reps=1, int gridsize = 100,std::string funname = "standard"){
  // OdeTraj is the one with low resolution
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,2) - OdeTraj.cols(1,2);
    //f_cur_centered(0,1)=0;
    //f_cur_centered(0,0)=0;
    //simulate a new trajectory
    List v = Traj_sim(state,OdeTraj,FTs);
    arma::mat v_traj = as<mat>(v[0]).cols(1,2) -  OdeTraj.cols(1,2);
    if(v_traj.has_nan()){
      Rcout<<v_traj<<endl;
      Rcout<<OdeTraj<<endl;
    }
    double u = R::runif(0,1);
  //  if(funname == "standard"){
      logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
 //   }else{
 //     logy = coal_loglik(init,f_cur,t_correct,lambda,gridsize) + log(u);
  //    }
    double theta = R::runif(0,2 * pi);

    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
    newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    int i = 0;
    while(newTraj.cols(1,2).min() <0 || coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize) <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){break;
        Rcout<<theta<<endl;}
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
      // newTraj.col(0) = f_cur.col(0);
      newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    }

    f_cur = newTraj;
   // Rcout<<coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize) <<endl;
    //Rcout<< logy - coal_loglik(init,newTraj,t_correct,lambda,gridsize)<<"1"<<endl;
  }
//Rcout<<coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize) <<endl;
  return newTraj;
}

//[[Rcpp::export()]]
arma::mat ESlice2(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                 List init, double t_correct, double lambda=10, int reps=1, int gridsize = 100,std::string funname = "standard",
                 bool volz = false, double beta = 0){
  // OdeTraj is the one with low resolution
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,2) - OdeTraj.cols(1,2);
    //f_cur_centered(0,1)=0;
    //f_cur_centered(0,0)=0;
    //simulate a new trajectory
    List v = Traj_sim(state,OdeTraj,FTs);
    arma::mat v_traj = as<mat>(v[0]).cols(1,2) -  OdeTraj.cols(1,2);
    if(v_traj.has_nan()){
      Rcout<<v_traj<<endl;
      Rcout<<OdeTraj<<endl;
    }
    double u = R::runif(0,1);
    //  if(funname == "standard"){
  //logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
   if(volz){
      logy = volz_loglik(init, LogTraj(f_cur), t_correct, beta, gridsize) + log(u);
    }else{
      logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
    }

    //   }else{
    //     logy = coal_loglik(init,f_cur,t_correct,lambda,gridsize) + log(u);
    //    }
    double theta = R::runif(0,2 * pi);

    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
    newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    int i = 0;
    double loglike;
    if(volz){
      loglike = volz_loglik(init, LogTraj(newTraj), t_correct, beta, gridsize);
    }else{
      loglike = coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize);
    }
    while(newTraj.cols(1,2).min() <0 || loglike <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){break;
        Rcout<<theta<<endl;}
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
      // newTraj.col(0) = f_cur.col(0);
      newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
      if(volz){
        loglike = volz_loglik(init, LogTraj(newTraj), t_correct, beta, gridsize);
      }else{
        loglike = coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize);
      }
    }

    f_cur = newTraj;
    // Rcout<<coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize) <<endl;
    //Rcout<< logy - coal_loglik(init,newTraj,t_correct,lambda,gridsize)<<"1"<<endl;
  }
  //Rcout<<coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize) <<endl;
  return newTraj;
}


