#include "basic_util.h"
#include "SIR_LNA.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

arma::mat SIRS_Fm_LNA_period(double X,double Y, double Z, double theta1,double theta2, double theta3,double alpha, double t,double period){
  arma::mat M;
  double th1 = theta1 * (1 + alpha * sin(t / period * 2 * pi));
  M << -th1 * Y << - th1 * X << theta3 <<arma::endr
    <<th1 * Y << th1 * X -theta2 <<0 <<arma::endr
    <<0 << theta2 << - theta3 <<endr;

  return M;
}

//[[Rcpp::export()]]
double DegenerateDet3(arma::mat M){
  return (M(1,1) - M(0,1)) * (M(2,2) - M(0,2)) - (M(1,2) - M(0,1)) * (M(1,2) - M(0,2));
}

//[[Rcpp::export()]]
arma::mat PseudoInverse(arma::mat M){
  double lambda1, lambda2;
  double delta = (M(1,1) + M(2,2) - M(0,1) - M(0,2)) * (M(1,1) + M(2,2) - M(0,1) - M(0,2)) -
    4 * ((M(1,1) - M(0,1)) * (M(2,2) - M(0,2)) - (M(1,2) - M(0,1)) * (M(1,2) - M(0,2)));
  if(delta<0.00000001){
    Rcout<<delta<<endl;
  }
  lambda1 = ((M(1,1) + M(2,2) - M(0,1) - M(0,2)) + sqrt(delta))/2;
  lambda2 = ((M(1,1) + M(2,2) - M(0,1) - M(0,2)) - sqrt(delta))/2;
  arma::vec eig1(3), eig2(3);
  eig1(0) = 1 - ((M(1,2) - M(0,1)) / (M(1,1) - M(0,1) - lambda1));
  eig1(1) = (M(1,2) - M(0,1)) / (M(1,1) - M(0,1) - lambda1);
  eig1(2) = -1;
  eig1 /= arma::norm(eig1);
  eig2(0) = 1 - ((M(1,2) - M(0,1)) / (M(1,1) - M(0,1) - lambda2));
  eig2(1) = (M(1,2) - M(0,1)) / (M(1,1) - M(0,1) - lambda2);
  eig2(2) = -1;
  eig2 /= arma::norm(eig2);
  arma::mat invs(3,3);
  invs.col(0) = eig1;
  invs.col(1) = eig2;
  invs.col(2) = sqrt(3) / 3 * arma::ones(3,1);
  arma::mat invmid = diagmat(arma::zeros(3));
  invmid(0,0) = 1 / lambda1;
  invmid(1,1) = 1 / lambda2;
  return invs * invmid * invs.t();
}
//[[Rcpp::export()]]
List SIRS_IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2,double theta3,double alpha,double period){
  arma::mat Sigma,F(3,3);
  Sigma.zeros(3,3);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(3);
  A << -1 << 1 << 0 <<arma::endr<<0<<-1<<1 <<arma::endr
    << 1 <<0 << -1<<arma::endr;
  arma::mat H,F0,Xinv;
  F0 = arma::diagmat(ones(3));
  double t;
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    double th1 = theta1 * (1 + alpha * sin(t / period * 2 * pi));
    h(0)= th1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = theta2 * Traj_par(i,2);
    h(2) = theta3 * Traj_par(i,3);
    H = diagmat(h);
    F = SIRS_Fm_LNA_period(Traj_par(i,1),Traj_par(i,2), Traj_par(i,3), theta1, theta2, theta3,alpha,t,period);
    F0 = F0 + (F * F0) * dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  A.t() * H * A ) * dt;
  }
  List Res;
  Res["expF"] = F0;
  Res["Simga"] = Sigma;
  if(Sigma.has_nan()){
    Rcout<<Traj_par<<endl;
  }
  return Res;
}

//[[Rcpp::export()]]
double log_like_trajSIRS(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                     int gridsize,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows - 1;
  int p = SdeTraj.n_cols - 1;
  double loglik = 0;
  //  int id1,id0;
  arma::vec Xd1, Xd0;

  Xd0 = Xd1 = (SdeTraj.submat(0,1,0,p)-OdeTraj.submat(0,1,0,p)).t();
  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
    //  id0 = i * gridsize;
    //  id1 = (i+1) * gridsize - 1;

    arma::mat SigInv = PseudoInverse(Scube.slice(i) );
    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),p)-OdeTraj.submat(i+1,1,i+1,p)).t();

    if(SdeTraj(i+1,0) <= t_correct){
      arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
      // Rcout <<  ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0)) <<endl;
      loglik += -log(DegenerateDet3(Scube.slice(i))) * 0.5 - 0.5 * INexp(0,0);
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

//[[Rcpp::export()]]
arma::vec SIRS_ODE(arma::vec states, arma::vec param,double t, double period){
  /**
   * One step ODE function
   * states: length 3 vector of (S, I, R)
   * param: length 4 vector of (beta, gamma, mu, A)
   * t: double, time
   * period: double, infection period
   *
   * return a vector of (dS,dI,dR)
   */
  double dx, dy, dz;
  double th1 = param[0] * (1 + param[3] * sin(2 * pi * t / period));
//double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
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
List SIRS_KOM_Filter(arma::mat OdeTraj, arma::vec param,int gridsize, double period){
  int n = OdeTraj.n_rows;
  double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  int p = OdeTraj.n_cols - 1;
  arma::cube Acube(p,p,k);
  arma::cube Scube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = SIRS_IntSigma(Traj_part,dt,param[0],param[1],param[2],param[3],period);
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
arma::mat ODE2(arma::vec initial, arma::vec t, arma::vec param, double period){
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
    k1 = SIRS_ODE(X0,param,t[i-1], period);
    k2 = SIRS_ODE(X0 + k1 * dt / 2, param, t[i-1], period);
    k3 = SIRS_ODE(X0 + k2 * dt / 2, param, t[i-1], period);
    k4 = SIRS_ODE(X0 + k3 * dt / 2, param, t[i-1], period);
    X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
    OdeTraj.submat(i,1,i,p) = X1.t();
  }
  return OdeTraj;
}


//[[Rcpp::export()]]
List Traj_sim_SIRS(arma::vec initial, arma::mat OdeTraj, List Filter,double t_correct){
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
    arma::mat noise = mvrnormArma2(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
  // Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj(i+1,0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * PseudoInverse(Sig) * noise;
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(DegenerateDet3(Sig))*0.5 + l1(0,0);
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


//[[Rcpp::export()]]
List Traj_sim_SIRS_ez(arma::vec initial, arma::vec times, arma::vec param,
                 int gridsize,double t_correct, double period){
  /**
   * Directly simulate the trajectory given initial value and paraemter value
   * initial: initial state
   * times: vector of time points for simulating ODE
   * gridsize: size of grid for LNA
   * t_correct: likelhooding cutting point for inference
   *
   * return a list contains the simulated trajectory and the likelihood
   *
   */
  int p = initial.n_elem;
  int k = times.n_rows / gridsize;
  arma::mat OdeTraj_thin = ODE2(initial,times,param, period);
  arma::mat OdeTraj(k+1,p+1);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
  }
  List Filter = SIRS_KOM_Filter(OdeTraj_thin,param,gridsize,period);

  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  k = OdeTraj.n_rows-1;

  arma::vec X0,X1 = initial,eta0(p),eta1 = initial;

  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = initial.t();
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
    eta1 = OdeTraj.submat(i+1, 1, i+1, p).t();
    arma::mat noise = mvrnormArma2(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj((i+1), 0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * PseudoInverse(Sig) * noise;
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(DegenerateDet3(Sig)) * 0.5 + l1(0,0);
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


//[[Rcpp::export()]]
arma::mat ESlice_SIRS(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                 List init, arma::vec params4, double t_correct, double lambda=10,
                 int reps=1, int gridsize = 100, bool volz = true){
  /**
   * Eliptical slice sampler for updating the latent trajectory
   *
   * params4: a length 4 vector of all parameters related to the betat, (beta, A, period, N)
   *
   */

  // OdeTraj is the one with low resolution
  int p = f_cur.n_cols - 1;
  arma::vec betat = betaDyn(params4[0], params4[1], f_cur.col(0), params4[2]) * params4[3];
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,p) - OdeTraj.cols(1,p);
    //f_cur_centered(0,1)=0;
    //f_cur_centered(0,0)=0;
    //simulate a new trajectory

    List v = Traj_sim_SIRS(state,OdeTraj,FTs,t_correct);

    arma::mat v_traj = as<mat>(v[0]).cols(1,p) -  OdeTraj.cols(1,p);
    if(v_traj.has_nan()){
      Rcout<<v_traj<<endl;
      Rcout<<OdeTraj<<endl;
    }
    double u = R::runif(0,1);
    //  if(funname == "standard"){
   // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
   if(volz){
     logy = volz_loglik_nh(init, LogTraj(f_cur), betat, t_correct, gridsize) + log(u);
     if(logy < -300){
       Rcout << log(u) <<endl;
     }
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
    newTraj.cols(1,p) = f_prime + OdeTraj.cols(1,p);
    int i = 0;
    double loglike;
    if(volz){
      loglike = volz_loglik_nh(init, LogTraj(newTraj),betat,t_correct, gridsize);
    }else{
      loglike = coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize);
    }
    while(newTraj.cols(1,p).min() <0 || loglike <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){
        Rcout<<theta<<endl;
        break;
        }
      if(theta < 0){
        theta_min = theta;
      }else{
        theta_max = theta;
      }
      theta = R::runif(theta_min,theta_max);
      f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta);
      // newTraj.col(0) = f_cur.col(0);

      newTraj.cols(1,p) = f_prime + OdeTraj.cols(1,p);

      if(volz){
        loglike = volz_loglik_nh(init, LogTraj(newTraj), betat, t_correct, gridsize);
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
