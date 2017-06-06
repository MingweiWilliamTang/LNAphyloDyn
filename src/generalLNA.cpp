#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

using namespace Rcpp;
using namespace arma;

//' @export Foo
class Foo{

  private:

    arma::mat A_pre; // pre reaction matrix
    arma::mat A_post; // post reaction matrix
    arma::mat A; // effect matrix in chemical reaction

  public:
    arma::vec param; // a vector of the parameters in reactions
    arma::ivec x_i; // a vector of other source of information
    arma::vec x_r; // a vection of other source of information
    int m; // number of reactions
    int p; // number of reactant
    double N;

    Foo(arma::mat A_pre_, arma::mat A_post_, arma::vec param_, arma::ivec x_i_, arma::vec x_r_){
      A_pre = A_pre_;
      A_post = A_post_;
      A = A_post - A_pre;
      param = param_;
      x_i = x_i_;
      x_r = x_r_;
      p = A_pre.n_cols;
      m = A_pre.n_rows;
      N = x_r[0];
    }

    // get functions
    arma::mat get_A_pre(){return A_pre;}
    arma::mat get_A_post() {return A_post;}
    arma::mat get_A(){return A;}
    arma::vec get_param() {return param;}

    // set functions
    void set_param(arma::vec new_param) {param = new_param;}
    void set_x_i(arma::ivec new_x_i ) {x_i = new_x_i;}
    void set_x_r(arma::vec new_x_r) {x_r = new_x_r;}
    //void set_A_pre(arma::mat new_A_pre) {A_pre = new_A_pre;}
    //void set_A_post(arma::mat new_A_post) {A_post = new_A_post;}


    arma::vec transform_param(double t){
      arma::vec param_new(m);
      param_new = param;
      param_new(0) = param(0) * param(1) / N;
      return param_new;
    }
/*
    arma::vec rate_h(arma::vec states, double t){
      arma::vec h(m);
      arma::vec param_new = transform_param(t);
      for(int i = 0; i < m; i ++){
        double s = param_new[i];
        for(int j = 0; j < p; j ++){
          if(A_pre[i,j] != 0){
          s *= (A_pre[i,j] * (A_pre[i,j] - 1) / 2) * arma::pow(states[j],A[i,j]);
          }
        }
        h(i) = s;
      }
      return h;
    }
*/

    arma::vec rate_h(arma::vec states, double t){
      arma::vec h(m);
      arma::vec param_new = transform_param(t);
      h(0) = param_new(0) * states(0) * states(1);
      h(1) = param_new(1) * states(1);
      h(2) = param_new(2) * N;
      return h;
    }
    arma::mat dh(arma::vec states, double t){
      arma::vec param_new = transform_param(t);
      arma::mat res;
      res << param_new[0] * states[1] << param_new[0] * states[0] << arma::endr
          << 0 << param_new[1] << arma::endr << 0.0 << 0.0 <<arma::endr;
      return res;
    }

    arma::vec F_vec(arma::vec states, double t){
     return A.t() * dh(states, t);
    }

    arma::vec ODE_ones(arma::vec states,double t){
      return A.t() * rate_h(states,t);
    }

    arma::mat ODE_rk45(arma::vec initial, arma::vec t){

      // XPtr<ODEfuncPtr> SIR_ODEfun = putFunPtrInXPtr(funname);
      //double N = initial[0] + initial[1] + initial[2];

      int n = t.n_rows;// number of time points on the grid

      double dt = t[1] - t[0];
      arma::mat OdeTraj(n,p + 1);
      OdeTraj.col(0) = t;
      OdeTraj.submat(0,1,0,p) = initial.t();
      arma::vec X0 = initial, k1=initial, k2=initial, k3=initial, k4=initial,X1=initial;
      for(int i = 1; i < n; i++){
        X0 = X1;
        k1 = ODE_ones(X0,t[i-1]);
        k2 = ODE_ones(X0 + k1 * dt / 2, t[i-1]);
        k3 = ODE_ones(X0 + k2 * dt / 2, t[i-1]);
        k4 = ODE_ones(X0 + k3 * dt / 2, t[i-1]);
        X1 = X0 + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
        OdeTraj.submat(i,1,i,p) = X1.t();
      }
      return OdeTraj;
    }

    List IntSigmaF(arma::mat Traj_par){

      arma::mat Sigma,F(2,2),F0(2,2);
      Sigma.zeros(2,2);
      int k = Traj_par.n_rows;
      arma::mat A;
      arma::vec h(2);
      arma::mat H,Xinv;
      F0 << 1 <<0 <<arma::endr << 0 << 1 <<endr;
      double t;
      double dt = Traj_par(1,0) - Traj_par(0,0);
      for(int i = 0; i < k; i++){
        t = Traj_par(i,0);
        arma::vec state(2);
        state(0) = Traj_par(i,1);
        state(1) = Traj_par(i,2);
        h = rate_h(state,t);
        H = diagmat(h);
        F = F_vec(state,t);
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

    List KOM_Filter_MX(arma::mat OdeTraj, int gridsize){
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
        List tempres = IntSigmaF(Traj_part);
        Acube.slice(i) = as<arma::mat>(tempres[0]);
        Scube.slice(i) = as<arma::mat>(tempres[1]);
      }
      List Res;
      Res["A"] = Acube;
      Res["Sigma"] = Scube;
      return Res;
    }

    double LNA_log_like_traj(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                                 int gridsize,double t_correct){
      arma::cube Acube = as<arma::cube>(Filter[0]);
      arma::cube Scube = as<arma::cube>(Filter[1]);

      int k = SdeTraj.n_rows - 1;
      // int p = SdeTraj.n_cols - 1;
      int p = 2;
      double loglik = 0;
      //  int id1,id0;
      arma::vec Xd1;
      arma::vec Xd0;

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
      }

      return loglik;
    }


    List LNA_Traj_sim(arma::mat OdeTraj, List Filter,double t_correct){
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

    List LNA_Traj_sim_ez(arma::vec initial, arma::vec times,
                      int gridsize,double t_correct){
      //int p = initial.n_elem;
      int p = 2;
      int k = times.n_rows / gridsize;
      arma::mat OdeTraj_thin = ODE_rk45(initial,times);
      arma::mat OdeTraj(k+1,p+1);
      for(int i = 0; i< k + 1; i++){
        OdeTraj.submat(i,0,i,p) = OdeTraj_thin.submat(i*gridsize,0,i*gridsize,p);
      }

      List Filter = KOM_Filter_MX(OdeTraj_thin, gridsize);

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


    arma::mat LNA_ESlice(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
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
        List v = LNA_Traj_sim(OdeTraj,FTs,t_correct);
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

        double theta_min = theta - 2 * pi;
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

};




RCPP_MODULE(mod_Foo){
  using namespace Rcpp;
  class_<Foo>( "Foo" )
  .constructor<arma::mat, arma::mat, arma::vec, arma::ivec, arma::vec>()
  .field( "x_i", &Foo::x_i)
  .field( "x_r", &Foo::x_r)
  .field( "param", &Foo::param)
  .field( "N", &Foo::N)
  .field( "m", &Foo::m)
  .field( "p", &Foo::p)
  .property( "param", &Foo::get_param, &Foo::set_param )
  .property( "A_pre", &Foo::get_A_pre)
  .property( "param", &Foo::get_A_post)
  .property( "param", &Foo::get_param)
  .method("get_param", &Foo::get_param)
  .method("get_A",&Foo::get_A)
  .method("IntSigmaF", &Foo::IntSigmaF)
  .method("KOM_Filter_MX",&Foo::KOM_Filter_MX)
  .method("set_param",&Foo::set_param)
  .method("LNA_Traj_sim", &Foo::LNA_Traj_sim)
  .method("LNA_log_like_traj", &Foo::LNA_log_like_traj)
  .method("set_x_i",&Foo::set_x_i)
  .method("set_x_r",&Foo::set_x_r)
  .method("get_A_pre",&Foo::get_A_pre)
  .method("get_A_post",&Foo::get_A_post)
  .method( "transform_param", &Foo::transform_param)
  .method( "rate_h", &Foo::rate_h)
  .method( "dh", &Foo::dh)
  .method( "F_vec", &Foo::F_vec)
  .method( "ODE_ones", &Foo::ODE_ones)
  .method( "ODE_rk45", &Foo::ODE_rk45)
  .method("LNA_Traj_sim_ez",&Foo::LNA_Traj_sim_ez)
  .method("LNA_ESlice",&Foo::LNA_ESlice);
  }
