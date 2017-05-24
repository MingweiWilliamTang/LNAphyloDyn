#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

//[[Rcpp::export()]]
arma::vec SIR_BD_period_ODE_one(arma::vec states, double N, arma::vec param, double t, double period){
  double dx, dy, dz;
  double th1 = param[0] * (1 + param[3] * sin(2 * pi * t / period));
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  dx = - th1 * states[0] * states[1] + param[2] * N - param[2] * states[0];
  dy = th1 * states[0] * states[1] - param[1] * states[1] -  param[2] * states[1];
  dz = param[1] * states[1] - param[2] * states[2];
  arma::vec res(3);
  res(0) = dx;
  res(1) = dy;
  res(2) = dz;
  return res;
}

//[[Rcpp::export()]]
arma::vec SIR_BD_period_ODE_one2(arma::vec states, double N, arma::vec param, double t, double period){
  double dx, dy;
  double th1 = param[0] * (1 + param[3] * sin(2 * pi * t / period));
  //double th1 = exp(param[0] + param[3] * sin(2 * pi * t / 40.0));
  dx = - th1 * states[0] * states[1] + param[2] * N - param[2] * states[0];
  dy = th1 * states[0] * states[1] - param[1] * states[1] -  param[2] * states[1];
  arma::vec res(2);
  res(0) = dx;
  res(1) = dy;
  return res;
}


//[[Rcpp::export()]]
arma::mat SIR_BD_Fm_LNA(double X,double Y, double theta1,double theta2, double mu,double alpha,
                        double t, double period){
  arma::mat M;
  double th1 = theta1 * (1 + alpha * sin(t / period * 2 * pi));
  M << - th1 * Y - mu<< - th1 * X <<arma::endr
    <<th1 * Y << th1 * X -theta2 - mu <<arma::endr;
  return M;
}

//[[Rcpp::export()]]
List SIR_BD_IntSigma(arma::mat Traj_par,double dt,double theta1,double theta2,double theta3,double alpha,
                     double N, double period){
  arma::mat Sigma,F(2,2),F0(2,2);
  Sigma.zeros(2,2);
  int k = Traj_par.n_rows;
  arma::mat A;
  arma::vec h(4);
  A << -1 << 1  <<arma::endr<<0<< -1 << arma::endr
    << 1 <<0  <<arma::endr << 0 << -1 << arma::endr;
  arma::mat H,Xinv;
  F0 << 1 <<0 <<arma::endr << 0 << 1 <<endr;
  double t;
  for(int i = 0; i < k; i++){
    t = Traj_par(i,0);
    double th1 = theta1 * (1 + alpha * sin(t / period * 2 * pi));
    h(0)= th1 * Traj_par(i,1) * Traj_par(i,2);
    h(1) = theta2 * Traj_par(i,2);
    h(2) = theta3 * (N - Traj_par(i,1));
    h(3) = theta3 * Traj_par(i,2);
    H = diagmat(h);
    F = SIR_BD_Fm_LNA(Traj_par(i,1),Traj_par(i,2), theta1, theta2, theta3,alpha,t,period);
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
arma::mat SIR_BD_ODE(arma::vec initial, arma::vec t, arma::vec param,
                     double N, double period){
  // XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
  //double N = initial[0] + initial[1] + initial[2];
  int n = t.n_rows;
  int p = initial.size();
  double dt = t[1] - t[0];
  arma::mat OdeTraj(n,p + 1);
  OdeTraj.col(0) = t;
  OdeTraj.submat(0,1,0,p) = initial.t();
  arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
  if(p == 2){
    for(int i = 1; i < n; i++){
      X0 = X1;
      k1 = SIR_BD_period_ODE_one2(X0,N,param,t[i-1], period);
      k2 = SIR_BD_period_ODE_one2(X0 + k1 * dt / 2, N,param, t[i-1], period);
      k3 = SIR_BD_period_ODE_one2(X0 + k2 * dt / 2, N,param, t[i-1], period);
      k4 = SIR_BD_period_ODE_one2(X0 + k3 * dt / 2, N,param, t[i-1], period);
      X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
      OdeTraj.submat(i,1,i,p) = X1.t();
    }
  }else{
      for(int i = 1; i < n; i++){
        X0 = X1;
        k1 = SIR_BD_period_ODE_one(X0,N,param,t[i-1],period);
        k2 = SIR_BD_period_ODE_one(X0 + k1 * dt / 2, N,param, t[i-1],period);
        k3 = SIR_BD_period_ODE_one(X0 + k2 * dt / 2, N,param, t[i-1],period);
        k4 = SIR_BD_period_ODE_one(X0 + k3 * dt / 2, N,param, t[i-1],period);
        X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
        OdeTraj.submat(i,1,i,p) = X1.t();
    }
  }
  return OdeTraj;
}



//[[Rcpp::export()]]
arma::mat SIR_BD_SDE(arma::vec init, double N, arma::vec param, arma::vec t, double period){
  int n = t.size();
  int p = init.size();
  arma::mat H;
  arma::mat Traj(n,p + 1);
  Traj.col(0) = t;
  Traj.submat(0,1,0,p) = init.t();
  arma::mat A;
  A << -1 << 1 << 0 << endr
    << 0 << -1 << 1 << endr
    << 1 <<  0 << 0 << endr
    << 0 << -1 << 0 << endr
    << 0 << 0 << -1 << endr;
  double dt  = t(1) - t(0);
  double beta = param(0);
  double gamma = param(1);
  double mu = param(2);
  double alpha = param(3);
  //double alpha = param(5);
  for( int i = 0; i < n-1; i ++){
    // double S_new, I_new, R_new;
     double betat = beta * (1 + alpha * (sin(2 * pi * t(i) / period)));
    // double betat2 = beta2 * (1 + alpha * (sin(2 * pi * t(i) / period)));
    arma::mat noises = randn(5,1) * sqrt(dt);
    arma::vec h(5);
    h(0) = betat * Traj(i,1) * Traj(i,2);
    h(1) = gamma * Traj(i,2);
    h(2) = (N - Traj(i,1)) * mu;
    h(3) = mu * Traj(i,2);
    h(4) = mu * Traj(i,3);
    //Rcout<<h<<endl;
    arma::mat mean = A.t() * h * dt;
    arma::mat change = mean + A.t() * arma::diagmat(arma::sqrt(h)) * noises;
    Traj.submat(i+1,1,i+1,p) = Traj.submat(i,1,i,p) + change.t();
  }
  return Traj;
}



//[[Rcpp::export()]]
arma::mat SIR_BD_period_SDE(arma::vec init, double N, arma::vec param, arma::vec t, double period ){
  int n = t.size();
  int p = init.size();
  arma::mat H;
  arma::mat Traj(n,p + 1);
  Traj.col(0) = t;
  Traj.submat(0,1,0,p) = init.t();
  double dt  = t(1) - t(0);
  double gamma = param(1);
  double alpha = param(2);
  double mu = param(3);
  double F = param(4);
  double beta = param(0) * (mu + gamma);
  for( int i = 0; i < n-1; i ++){
    // double S_new, I_new, R_new;
    double betat = beta * (1 + alpha * (sin(2 * pi * t(i) / period)));
    arma::mat noises = randn(2,1);
    Traj(i+1,1) = Traj(i,1) + (mu * N - mu*Traj(i,1) - (betat*Traj(i,1)*Traj(i,2)/N)) * dt +
      F*(betat*Traj(i,1)*Traj(i,2)/N) * sqrt(dt) * noises(0,0);
    Traj(i+1,2) = Traj(i,2) + (betat*Traj(i,1)*Traj(i,2)/N - gamma*Traj(i,2) - mu*Traj(i,2)) * dt +
      F*(betat*Traj(i,1)*Traj(i,2)/N) * sqrt(dt) * noises(1,0);
    Traj(i+1,3) = Traj(i,3) + (gamma * Traj(i,2) - mu*Traj(i,3)) * dt;
  }
  return Traj;
}



//[[Rcpp::export()]]
List SIR_BD_KOM_Filter(arma::mat OdeTraj, arma::vec param,int gridsize, double N, double period){
  int n = OdeTraj.n_rows;
  double dt = (OdeTraj(1,0) - OdeTraj(0,0));
  int k = (n-1) / gridsize;
  //int p = OdeTraj.n_cols - 1;
  int p =2;
  arma::cube Acube(p,p,k);
  arma::cube Scube(p,p,k);
  arma::mat Traj_part;
  for(int i=0;i<k;i++){
    Traj_part = OdeTraj.submat(i*gridsize,0,(i+1)*gridsize-1,p);
    List tempres = SIR_BD_IntSigma(Traj_part,dt,param[0],param[1],param[2],param[3], N,period);
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
double log_like_trajSIR_BD(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                           int gridsize,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);

  int k = SdeTraj.n_rows - 1;
 // int p = SdeTraj.n_cols - 1;
  int p = 2;
  double loglik = 0;
  //  int id1,id0;
  arma:vec Xd1, Xd0;

  Xd0 = Xd1 = (SdeTraj.submat(0,1,0,p)-OdeTraj.submat(0,1,0,p)).t();
  // Xd0 = Xd1 = (SdeTraj.submat(0,1,0,3) - OdeTraj.submat(0,1,0,3)).t();
  for(int i = 0; i < k; i++){
    Xd0 = Xd1;
    //  id0 = i * gridsize;
    //  id1 = (i+1) * gridsize - 1;

    arma::mat SigInv = inv2(Scube.slice(i) );
    arma::mat A = Acube.slice(i);

    Xd1 = (SdeTraj.submat((i+1),1,(i+1),p)-OdeTraj.submat(i+1,1,i+1,p)).t();

    if(SdeTraj(i+1,0) <= t_correct && i != 0){
      arma::mat INexp = ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0));
      // Rcout <<  ((Xd1 - A * Xd0).t() * SigInv * (Xd1 - A * Xd0)) <<endl;
      loglik += -log(arma::det(Scube.slice(i))) / 2 - 0.5 * INexp(0,0) ;
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
List Traj_sim_SIR_BD(arma::vec initial, arma::mat OdeTraj, List Filter,double t_correct){
  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  int k = OdeTraj.n_rows - 1;
  // int p = OdeTraj.n_cols - 1;
  int p = 2;
  arma::vec X0,X1 = initial.subvec(0,p-1),eta0(p),eta1 = initial.subvec(0,p-1);
  double loglike = 0;
  arma::mat LNA_traj(k+1,p+1);
  LNA_traj(0,0) = 0;

  LNA_traj.submat(0,1,0,p) = initial.subvec(0,p-1).t();
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
List Traj_sim_SIR_BD_ez(arma::vec initial, arma::vec times, arma::vec param,
                        int gridsize, double N,double t_correct,
                        double period){
  //int p = initial.n_elem;
  int p = 2;
  int k = times.n_rows / gridsize;
  arma::mat OdeTraj_thin = SIR_BD_ODE(initial,times,param,N, period);
  arma::mat OdeTraj(k+1,p+1);
  for(int i = 0; i< k + 1; i++){
    OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
  }
  List Filter = SIR_BD_KOM_Filter(OdeTraj_thin,param,gridsize,N,period);

  arma::cube Acube = as<arma::cube>(Filter[0]);
  arma::cube Scube = as<arma::cube>(Filter[1]);
  k = OdeTraj.n_rows-1;
  //Rcout << Scube.slice((k-1)) <<endl;
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
arma::mat ESlice_SIR_BD(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
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
    List v = Traj_sim_SIR_BD(state,OdeTraj,FTs,t_correct);
    arma::mat v_traj = as<mat>(v[0]).cols(1,p) -  OdeTraj.cols(1,p);
    if(v_traj.has_nan()){
      Rcout<<v_traj<<endl;
      Rcout<<OdeTraj<<endl;
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
        Rcout <<"vvvvv"<<endl;
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


//[[Rcpp::export()]]
double StateSpace_Like(arma::vec state,double beta, double N, arma::vec w, arma::vec C, arma::vec y){
  arma::vec ll = -2 * (state(0) / state(1)) * beta * N * (w % C) + \
    (log(state(0)) - log(state(1)) + log(beta * N)) * y;
  return exp(sum(ll));
}
//[[Rcpp::export()]]
double StateSpace_Period_Like2(arma::vec state,double beta, double N, arma::vec w, arma::vec C, arma::vec y,double t,double A, double period){
  double betat = beta * (1 + A * sin(2*pi * t));
  arma::vec ll = -2 * (state(0) / state(1)) * betat * N * (w % C) + \
    (log(state(0)) - log(state(1)) + log(betat * N)) * y;
  return exp(sum(ll));
}

//[[Rcpp::export()]]
arma::vec SampleWithReplace(arma::vec prob, int N){
  int d = prob.n_elem;
  arma::vec res(N);
  for(int i = 0; i<N; i++){
    double q = R::runif(0,1);
    double cumsum = 0;
    for(int j = 0; j < d; j++){
      cumsum += prob(j);
      if(q < cumsum){
        res(i) = j;
        break;
      }
    }
  }
  return res;
}



//[[Rcpp::export()]]
List SIR_BD_SMC(arma::vec params,double N, List init, int D, arma::vec TimeGrid, int OdeSize, double t_correct,arma::vec prior){
  // D is the number of the paraticles
  // init is the list for measuring colescent likelihood
  int n0 = 0; // the index of the last time point
//  double alpha = params(0);
  double beta0 = params(0);
  double gamma = params(1);
  double mu = params(2);
  double A = params(3);
  double period = params(4);
  arma::vec w = as<vec>(init[3]);
  arma::vec C = as<vec>(init[2]);
  arma::vec gridrep = as<vec>(init[6]);
  arma::vec y = as<vec>(init[4]);
  arma::vec t_all = as<vec>(init[0]);
  while(TimeGrid(n0) <= t_all.max()){
    n0 ++;
  }
  // Rcout << n0 << endl;
  arma::vec grid_idx = as<vec>(init[12]);
  // Rcout << grid_idx.n_elem <<endl;
  double dt = (TimeGrid(1) - TimeGrid(0)) / OdeSize;
  arma::vec wi(D);
  arma::vec W = arma::ones(D)/D;
  arma::vec State(2);
  arma::cube Trajs(n0,3,D);
  arma::mat TrajTemp(D,2);// particle trajs
  arma::mat OdeTrajTemp(D,2);
  arma::mat OdeTraj(D,2);
  arma::vec py(n0-1);
  arma::vec Pid;
  arma::vec ts(OdeSize + 1);
  double t0 = TimeGrid(1);
  for(int i=1; i < n0; i++){
    // resample
    if(i > 1){
    Pid = SampleWithReplace(W,D);
      //Pid = arma::ones(D)/D;
    for(int k = 0; k < OdeSize; k++){
      ts(k) = t0 + k * dt;
    }
      ts(OdeSize) = ts(OdeSize - 1) + dt;
      t0 = ts(OdeSize);
    }
      int id1 = int(grid_idx(n0-i-1));
      int id2 = int(grid_idx(n0-i))-1;
      //Rcout<<id1 <<"\t"<<id2 <<endl;
    for(int d = 0; d < D; d ++){
      if(i == 1){
        double alpha =  exp(R::runif(prior(0),prior(1)));
        State(0) = N * alpha / (1 + alpha);
        State(1) = N - State(0);
       // Rcout << State<< endl;
        TrajTemp.submat(d,0,d,1) = State.t();
        OdeTrajTemp.submat(d,0,d,1) = State.t();
        // wi(d) = StateSpace_Like(State,beta0,N,w.subvec(id1,id2),
          //  C.subvec(id1,id2),y.subvec(id1,id2));
        wi(d) = StateSpace_Period_Like2(State,beta0,N,w.subvec(id1,id2),
            C.subvec(id1,id2),y.subvec(id1,id2),t0,A,period);

         Trajs.slice(d).submat(i,0,i,0) = t0;
      }else{
        //Rcout<<i<<endl;

        Trajs.slice(d).submat(i,0,i,0) = t0;
        Trajs.slice(d).submat(i-1,1,i-1,2) =  TrajTemp.submat(int(Pid(d)),0,int(Pid(d)),1);
        OdeTraj.submat(d,0,d,1) = OdeTrajTemp.submat(int(Pid(d)),0,int(Pid(d)),1);
        arma::mat OneStepODE = SIR_BD_ODE(OdeTraj.submat(d,0,d,1).t(),
                                            ts,params.subvec(0,3), N,period);
        List tempres = SIR_BD_IntSigma(OneStepODE,dt,beta0,gamma,mu,A, N,period);
        OdeTrajTemp.submat(d,0,d,1) = OneStepODE.submat(OdeSize,1,OdeSize,2);
        arma::mat F = as<arma::mat>(tempres[0]);
        arma::mat Sigma = as<arma::mat>(tempres[1]);
        arma::mat noise = mvrnormArma(2,Sigma);
        State = F * (Trajs.slice(d).submat(i-1,1,i-1,2) - OdeTraj.submat(d,0,d,1)).t() +
        OdeTrajTemp.submat(d,0,d,1).t() + noise;
        TrajTemp.submat(d,0,d,1) = State.t();
        wi(d) = wi(d) = StateSpace_Period_Like2(State,beta0,N,w.subvec(id1,id2),
           C.subvec(id1,id2),y.subvec(id1,id2),t0,A,period);
        if(i == n0 - 1){
          Trajs.slice(d).submat(n0-1,1,n0-1,2) = TrajTemp.submat(d,0,d,1);
        }
      }
      }
    py(i-1) = log(sum(wi) / D);
    //if(py.has_nan()){
    //  Rcout<< i<<endl;
  //  }
 if(log(sum(wi)) < -100){
   py = arma::ones(D) * (-100000000);
   break;
 }
  }

  List res;
  res["loglike"] = sum(py);
  res["Traj"] = Trajs;
  arma::mat MeanTraj(n0-1,3);

  MeanTraj.col(0) = Trajs.slice(0).col(0).subvec(1,n0-1);
  for(int i = 0; i < n0 - 1; i++){
    double I = 0;
    double S = 0;
    for(int j = 0; j < D; j++){
      S += W(j) * Trajs(i+1,1,j);
      I += W(j) * Trajs(i+1,2,j);
    }
    MeanTraj(i,1) = S;
    MeanTraj(i,2) = I;
  }
  res["MeanTraj"] = MeanTraj;
  return res;
}





//[[Rcpp::export()]]
List SIR_BD_PMCMC(List Init, double N,int D, arma::vec TimeGrid, int OdeSize, double t_correct, int niter,
                  double period,double mu,double A,
                  arma::vec priorAlpha, arma::vec priorGamma,double pR0,double pgamma){
// initialization
double R0 = R::runif(1,10);
double gamma = exp(R::rnorm(priorGamma(0),priorGamma(1)));
double logGamma = R::dnorm4(log(gamma),priorGamma(0),priorGamma(1),1);
Rcout << logGamma <<endl;
double beta = R0 * gamma / N;
arma::vec params(5);
params(0) = beta;
params(1) = gamma;
params(2) = mu;
params(3) = A;
params(4) = period;
List SMC_res = SIR_BD_SMC(params,N, Init, D,  TimeGrid, OdeSize, t_correct, priorAlpha);
double coal_lik = double(SMC_res["loglike"]);
arma::vec param_new = params;
arma::mat MCMC_store(niter,6);
int n0 = 0;
arma::vec t_all = as<vec>(Init[0]);
while(TimeGrid(n0) <= t_all.max()){
  n0 ++;
}

arma::cube TrajCube(n0-1,3,niter);
Rcout << niter <<endl;
for(int i=1; i<niter; i++){
  double R0_new = R0 + R::runif(-pR0,pR0);
  double gamma_new = gamma * exp(R::rnorm(pgamma,0.5));
  double beta_new = R0_new * gamma_new / N;
  param_new(0) = beta_new;
  param_new(1) = gamma_new;
  double logGamma_new = R::dnorm4(log(gamma_new),priorGamma(0),priorGamma(1),1);
  int AR1 = 0;
  if(R0_new < 1 || R0_new > 10) {
    MCMC_store.submat(i,0,i,3) = params.subvec(0,3).t();
    MCMC_store(i,4) = AR1;
    MCMC_store(i,5) = coal_lik;
    continue;
  }

    List SMC_new = SIR_BD_SMC(param_new,N, Init, D,  TimeGrid, OdeSize, t_correct, priorAlpha);
    double coal_lik_new = double(SMC_new["loglike"]);

    if(log(R::runif(0,1)) < logGamma_new - logGamma + coal_lik_new - coal_lik){

      AR1 = 1;
      logGamma = logGamma_new;
      gamma = gamma_new;
      R0 = R0_new;
      coal_lik = coal_lik_new;
      params(0) = beta_new;
      params(1) = gamma_new;
      SMC_res = SMC_new;
    }
    //Rcout << "ss" <<endl;
    MCMC_store.submat(i,0,i,3) = params.subvec(0,3).t();
    MCMC_store(i,4) = AR1;
    MCMC_store(i,5) = coal_lik;
    TrajCube.slice(i) = as<arma::mat>(SMC_res["MeanTraj"]);
    if(i % 100 == 0){
      Rcout << "finish "<< i <<" iterations" <<endl;
    }
  }
  List Res;
  Res["LatentTraj"] = TrajCube;
  Res["pars"] = MCMC_store;
  return Res;
}



/*

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
*/
