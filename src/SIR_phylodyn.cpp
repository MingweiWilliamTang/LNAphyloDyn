//@copyright Mingwei
//#include<RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
#include "basic_util.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
#include<iostream>
using namespace Rcpp;
using namespace arma;
# define pi           3.1415926535897932384

arma::mat Fm_log_LNA(double X,double Y,double theta1,double theta2){
  arma::mat M;
  M << theta1 * exp(Y-X)/2 << -theta1*exp(Y)-theta1/2*exp(Y-X) <<arma::endr
    <<theta1*exp(X)-theta1*exp(X-Y)/2 << theta1*exp(X-Y)/2+theta2/2*exp(-Y)<<arma::endr;
  return M;
}

double coal_log_like_hetero(arma::vec ts, arma::vec traj, List gene){
  return 1.0;
}

//[[Rcpp::export()]]
double log_like_traj(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                     int gridsize,double t_correct = 90){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows-1;
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

    if(SdeTraj(i+1,0) <= t_correct){
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
List SIR_log_KOM_Filter2(arma::mat OdeTraj,double theta1,double theta2,int gridsize){
  int n = OdeTraj.n_rows;
  double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  arma::cube Acube(2,2,k);
  arma::cube Scube(2,2,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
//    Rcout<<i<<endl;
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,2);
    List tempres(IntSigma(Traj_part,dt,theta1,theta2));
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

//[[Rcpp::export()]]
arma::mat ODE(arma::vec initial, arma::vec t, arma::vec param){
  int n = t.n_rows;
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,3);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,2) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  for(int i = 1; i < n; i++){
    X0 = X1;
    k1 = SIR_ODE(X0[0],X0[1],param[0],param[1]);
    k2 = SIR_ODE(X0[0]+k1[0]*dt/2,X0[1]+k1[1]*dt/2,param[0],param[1]);
    k3 = SIR_ODE(X0[0]+k2[0]*dt/2,X0[1]+k2[1]*dt/2,param[0],param[1]);
    k4 = SIR_ODE(X0[0]+k3[0]*dt/2,X0[1]+k3[1]*dt/2,param[0],param[1]);
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
//    Rcout<<"test2"<<endl;
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
        int gridsize,double t_correct = 90){
  int k = times.n_rows / gridsize;
  arma::vec param(2);
  param(0) = theta1;
  param(1) = theta2;
  arma::mat OdeTraj_thin = ODE(initial,times,param);
  arma::mat OdeTraj(k+1,3);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,2) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,2);
  }
  List Filter = SIR_log_KOM_Filter2(OdeTraj_thin,theta1,theta2,gridsize);

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
    //   arma::mat SigInv = inv2(Scube.slice(i).submat(0,0,1,1));
    arma::mat A = Acube.slice(i);

    X0 = X1;

    eta0 = eta1;
    eta1 = OdeTraj.submat(i+1, 1, i+1, 2).t();
    //    Rcout<<"test2"<<endl;
    arma::mat noise = mvrnormArma(2,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
   if(OdeTraj((i+1), 0) <= t_correct){
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

/*
 * In put a SDE path trajectory of log-transformed SIR and the parameter value,
 return the loglikelihood
 */
//[[Rcpp::export()]]
double log_like_traj2(arma::mat SdeTraj,arma::vec times,arma::vec state,
                      double theta1,double theta2,int gridsize,double t_correct = 90){

// generate the ODE path, calculate mean and covariance in the transition probability
  int k = SdeTraj.n_rows-1;
  arma::vec param(2);
  param(0) = theta1;
  param(1) = theta2;

  arma::mat OdeTraj_thin = ODE(state,times,param);
  arma::mat OdeTraj(k+1,3);
  for(int i = 0; i< SdeTraj.n_rows; i++){
    OdeTraj.submat(i,0,i,2) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,2);
  }
  List Filter = SIR_log_KOM_Filter2(OdeTraj_thin,theta1,theta2,gridsize);
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
double coal_loglik(List init, arma::mat f1, double t_correct, double lambda, int gridsize = 1){
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
  while(tf(n0) < t_correct){
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
double coal_loglik_step(List init, arma::mat f1, double t_correct, double lambda, int gridsize = 1){
  int n0 = 0;
  while(f1(n0,0) < t_correct){
    n0 ++;
  }
  //  Rcout<<n0<<endl;
  arma::vec f2(n0);
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




/*
 * Input the current trajectory and the ODE trajectory as well as the
 *
 */
//[[Rcpp::export()]]
arma::mat ESlice(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                    List init, double t_correct, double lambda=10, int reps=1, int gridsize = 100){
  // OdeTraj is the one with low resolution
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,2) - OdeTraj.cols(1,2);

    //simulate a new trajectory
    List v = Traj_sim(state,OdeTraj,FTs);
    arma::mat v_traj = as<mat>(v[0]).cols(1,2) -  OdeTraj.cols(1,2);

    double u = R::runif(0,1);
    double logy = coal_loglik(init,f_cur,t_correct,lambda,gridsize) + log(u);

    double theta = R::runif(0,2 * pi);
    double theta_min = theta - 2*pi;
    double theta_max = theta;

    arma::mat f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
    newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    while(coal_loglik(init,newTraj,t_correct,lambda,gridsize) <= logy){
    // shrink the bracket

    if(theta < 0){
      theta_min = theta;
    }else{
      theta_max = theta;
    }
    theta = R::runif(theta_min,theta_max);
    f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
   // newTraj.col(0) = f_cur.col(0);
    newTraj.cols(1,2) = f_prime + OdeTraj.cols(1,2);
    f_cur = newTraj;
    }
  }
  return newTraj;
}

/*
List SampleS1(double s1,double s2,double state,double theta1,double theta2,
              double lambda,bool init = true){

  double s1_new;
  if(init == true){
    s1_new = s1 * exp(R::runif(-1,1));
  }else{
    s1_new = s1 + R::runif(-0.5,0.5) * 0.000002;
  }
  arma::vec param_new(2);
  param_new(0) = s1_new;
  param_new(1) = s1_new * s2;
  arma::mat Ode_Traj_thin_new = ODE(log(state),times,
                           param_new);
}

 */



