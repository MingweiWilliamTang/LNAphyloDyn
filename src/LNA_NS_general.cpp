#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

//[[Rcpp::export()]]
double betaf(double t, arma::vec param, arma::vec x_r, arma::ivec x_i){
/**
* x_r = (N,cht1, cht2, ...)
* theta = (R0,gamma, lambda, ch1,ch2,...)
* x_i = (nch,nparam)
*/
  double R0 = param[0];
  int i = 0;
  int nch = x_i[0];
  while(x_r[i + 1] <= t){
    R0 *= param[x_i[1] + i];
    if(i == nch - 1) break;
      i ++;
    }
    return R0 * param[1] / x_r[0];
}


//[[Rcpp::export()]]
arma::vec betafs(arma::vec ts, arma::vec param, arma::vec x_r, arma::ivec x_i){
  /**
   * x_r = (N,cht1, cht2, ...)
   * theta = (R0,gamma, lambda, ch1,ch2,...)
   * x_i = (nch,nparam)
   */
  double R0 = param[0];
  int i = 0;
  int nch = x_i[0];
  int m = ts.n_elem;
  arma::vec betas(m);
  for(int j = 0; j < m; j++){
    if(i < nch){
      while(x_r[i + 1] <= ts[j]){
        R0 *= param[x_i[1] + i];
        i ++;
        if(i == nch) break;
      }
    }
    betas(j) = R0 * param[1] / x_r[0];
  }
  return betas;
}



//[[Rcpp::export()]]
arma::vec ODE_general_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i){
  double dx, dy ;
  double th1 = betaf(t,param, x_r,x_i);
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  dx = - th1 * states[0] * states[1];
  dy = th1 * states[0] * states[1] - param[1] * states[1];
  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}


//[[Rcpp::export()]]
arma::mat general_F(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i){
  arma::mat M;
   double th1 = betaf(t,param, x_r,x_i);
  M << - th1 * states[1] << - th1 * states[0] <<arma::endr
    <<th1 * states[1] << th1 * states[0] -param[1] <<arma::endr;
  return M;
}


//[[Rcpp::export()]]
arma::mat General_ODE_rk45(arma::vec initial, arma::vec t, arma::vec param,
                      arma::vec x_r, arma::ivec x_i){
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
      k1 = ODE_general_one(X0,param,t[i-1], x_r,x_i);
      k2 = ODE_general_one(X0 + k1 * dt / 2,param, t[i-1], x_r,x_i);
      k3 = ODE_general_one(X0 + k2 * dt / 2,param, t[i-1], x_r,x_i);
      k4 = ODE_general_one(X0 + k3 * dt / 2,param, t[i-1], x_r,x_i);
      X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
      OdeTraj.submat(i,1,i,p) = X1.t();
    }
  return OdeTraj;
}



//[[Rcpp::export()]]
List General_IntSigmaF(arma::mat Traj_par,arma::vec param,
                     arma::vec x_r, arma::ivec x_i){
  arma::mat Sigma,F(2,2),F0(2,2);
  Sigma.zeros(2,2);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(2);
  A << -1 << 1  <<arma::endr<<0<< -1 << arma::endr;
  arma::mat H,Xinv;
  F0 << 1 <<0 <<arma::endr << 0 << 1 <<endr;
  double t;
  double dt = Traj_par(1,0) - Traj_par(0,0);
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    double th1 = betaf(t,param, x_r,x_i);
    arma::vec state(2);
    state(0) = Traj_par(i,1);
    state(1) = Traj_par(i,2);
    h(0)= th1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = param[1] * Traj_par(i,2);
    H = diagmat(h);
    F = general_F(state, param,t, x_r, x_i);
    F0 = F0 + (F * F0) * dt;
    Sigma = Sigma + (Sigma * F.t() + F * Sigma +  A.t() * H * A ) * dt;
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


//[[Rcpp::export()]]
List General_KOM_Filter(arma::mat OdeTraj, arma::vec param,int gridsize,arma::vec x_r, arma::ivec x_i){
  int n = OdeTraj.n_rows;
  //double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  //int p = OdeTraj.n_cols - 1;
  int p =2;
  arma::cube Acube(p,p,k);
  arma::cube Scube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = General_IntSigmaF(Traj_part,param,x_r,x_i);
    // Rcout<< tempres<<endl;
    Acube.slice(i) = as<arma::mat>(tempres[0]);
    Scube.slice(i) = as<arma::mat>(tempres[1]);
  }
  List Res;
  Res["A"] = Acube;
  Res["Sigma"] = Scube;
  return Res;
}



//[[Rcpp::export()]]
double log_like_traj_general(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                           int gridsize,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows - 1;
  // int p = SdeTraj.n_cols - 1;
  int p = 2;
  double loglik = 0;
  //  int id1,id0;
  arma::vec Xd1, Xd0;

  Xd1 = (SdeTraj.submat(0,1,0,p)-OdeTraj.submat(0,1,0,p)).t();
  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
    //  id0 = i * gridsize;
    //  id1 = (i+1) * gridsize - 1;

    arma::mat SigInv = inv2(Scube.slice(i) );
    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),p)-OdeTraj.submat(i+1,1,i+1,p)).t();

    if(SdeTraj(i+1,0) <= t_correct){
      arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
      // Rcout <<  ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0)) <<endl;
      loglik += -log(arma::det(Scube.slice(i)))/2.0 - 0.5 * INexp(0,0);
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
List Traj_sim_general(arma::mat OdeTraj, List Filter,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  int k = OdeTraj.n_rows - 1;
  // int p = OdeTraj.n_cols - 1;
  int p = 2;
  arma::vec X0(p),eta0(p),X1(p),eta1(p);
  X1(0) = OdeTraj(0,1);
  X1(1) = OdeTraj(0,2);
  eta1 = X1;
  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = X1.t();
  //Rcout<<"test1"<<endl;

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
      arma::mat l1 = (-0.5) * noise.t() * inv2(Sig) * noise;
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(arma::det(Sig))/2 + l1(0,0);
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
List Traj_sim_ezG(arma::vec initial, arma::vec times, arma::vec param,
                        int gridsize,arma::vec x_r,arma::ivec x_i,double t_correct){
  //int p = initial.n_elem;
  int p = 2;
  int k = times.n_rows / gridsize;
  arma::mat OdeTraj_thin = General_ODE_rk45(initial,times, param, x_r, x_i);
  arma::mat OdeTraj(k+1,p+1);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
  }

  List Filter = General_KOM_Filter(OdeTraj_thin, param, gridsize, x_r, x_i);

  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  k = OdeTraj.n_rows-1;
  //Rcout << Acube.slice(9) <<endl;
  arma::vec X0,X1 = initial.subvec(0,1),eta0(p),eta1 = initial.subvec(0,1);

  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = initial.subvec(0,1).t();
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
    arma::mat noise = mvrnormArma(p,Sig);

    X1 = A * (X0 - eta0) + eta1 + noise;
    //Rcout<<noise.submat(0,0,1,0)<<endl;
    if(OdeTraj((i+1), 0) <= t_correct){
      arma::mat l1 = (-0.5) * noise.t() * inv2(Sig) * noise;
      //     arma::mat l1  = (-0.5) * noise.t() * inv2(Sig + eye(3,3)*0.00001) * noise;
      //    loglike += -1.5 * log(det(Sig+eye(3,3)*0.00001)) + l1(0,0);
      loglike += -log(arma::det(Sig))/2 + l1(0,0);
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
arma::mat ESlice_general(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                        List init, arma::vec betaN, double t_correct, double lambda=10,
                        int reps=1, int gridsize = 100, bool volz = false){
  // OdeTraj is the one with low resolution
  int p = f_cur.n_cols - 1;
  arma::mat newTraj(f_cur.n_rows, f_cur.n_cols);
  double logy;
  for(int count = 0; count < reps; count ++){
    // centered the old trajectory without time grid
    arma::mat f_cur_centered = f_cur.cols(1,p) - OdeTraj.cols(1,p);
    //f_cur_centered(0,1)=0;
    //f_cur_centered(0,0)=0;
    //simulate a new trajectory
    List v = Traj_sim_general(OdeTraj,FTs,t_correct);
    arma::mat v_traj = as<mat>(v[0]).cols(1,p) -  OdeTraj.cols(1,p);
    if(v_traj.has_nan()){
      Rcout<<"dddd"<<endl;
    }
    double u = R::runif(0,1);
    //  if(funname == "standard"){
    // logy = coal_loglik(init,LogTraj(f_cur),t_correct,lambda,gridsize) + log(u);
    if(volz){
      logy = volz_loglik_nh(init, LogTraj(f_cur), betaN, t_correct, gridsize) + log(u);
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
      loglike = volz_loglik_nh(init, LogTraj(newTraj),betaN,t_correct, gridsize);
    }else{
      loglike = coal_loglik(init,LogTraj(newTraj),t_correct,lambda,gridsize);
    }
    while(newTraj.cols(1,p).min() <0 || loglike <= logy){
      // shrink the bracket
      i += 1;
      if(i>20){
        newTraj = f_cur;
        Rcout<<"theta = "<<theta<<endl;
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
        loglike = volz_loglik_nh(init, LogTraj(newTraj), betaN, t_correct, gridsize);
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

