#ifndef LNA_FUNCTIONAL_H_
#define LNA_FUNCTIONAL_H_

#include "basic_util.h"
#include "SIR_phylodyn.h"
#define pi 3.14159265358979323846264338327950288

typedef arma::vec (*parat)(double, arma::vec, arma::vec, arma::ivec);
typedef arma::vec (*ODE_fun)(arma::vec, arma::vec, double, arma::vec, arma::ivec, std::string, std::string);
typedef arma::mat (*F_fun)(arma::vec, arma::vec, std::string);
typedef arma::vec (*h_fun)(arma::vec, arma::vec, std::string);



arma::vec betaTs(arma::vec param, arma::vec times, arma::vec x_r, arma::ivec x_i);
arma::vec param_transform(double t, arma::vec param, arma::vec x_r, arma::ivec x_i);
XPtr<parat> transformPtr(std::string trans);
arma::vec ODE_SIR_one(arma::vec states, arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
            std::string transP , std::string transX);
arma::vec ODE_SEIR2_one(arma::vec states, arma::vec param, double t,
              arma::vec x_r, arma::ivec x_i,
              std::string transP,
              std::string transX);


arma::vec ODE_SEIR_one(arma::vec states, arma::vec param, double t,
             arma::vec x_r, arma::ivec x_i,
             std::string transP,
             std::string transX);


arma::vec ODE_SIRS_one(arma::vec states,arma::vec param, double t, arma::vec x_r, arma::ivec x_i,
             std::string transP, std::string transX);

arma::mat SIRS_F(arma::vec states,arma::vec thetas,std::string transX);

arma::mat SIR_F(arma::vec states,arma::vec thetas,std::string transX);

arma::mat SEIR2_F(arma::vec states, arma::vec thetas, std::string transX);

arma::vec SIR_h(arma::vec states,arma::vec thetas,std::string transX);


arma::vec SEIR_h(arma::vec states,arma::vec thetas,std::string transX);

XPtr<F_fun> F_funPtr(std::string model);

XPtr<h_fun> h_fun_Ptr(std::string model);

XPtr<ODE_fun> ODE_fun_Ptr(std::string model);


arma::mat ODE_rk45(arma::vec initial, arma::vec t, arma::vec param,
                   arma::vec x_r, arma::ivec x_i,
                   std::string transP,std::string model,
                   std::string transX);

List SigmaF(arma::mat Traj_par,arma::vec param,
            arma::vec x_r, arma::ivec x_i,
            std::string transP, std::string model,std::string transX);


List KF_param(arma::mat OdeTraj, arma::vec param,int gridsize,arma::vec x_r, arma::ivec x_i,
              std::string transP,
              std::string model,std::string transX);

List KF_param_chol(arma::mat OdeTraj, arma::vec param,int gridsize,arma::vec x_r, arma::ivec x_i,
                   std::string transP,
                   std::string model,std::string transX);

double log_like_traj_general2(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                              int gridsize,double t_correct);

double log_like_traj_general_ez(arma::mat SdeTraj, double t_correct, arma::vec initial, arma::vec t, arma::vec param,
                                arma::vec x_r, arma::ivec x_i,
                                std::string transP,std::string model,
                                std::string transX);


double log_like_traj_general_adjust(arma::mat SdeTraj,arma::mat OdeTraj, List Filter_NC,
                                    int gridsize,double t_correct);



List Traj_sim_general3(arma::mat OdeTraj, List Filter,double t_correct);


List Traj_sim_general_noncentral(arma::mat OdeTraj, List Filter_NC,double t_correct);

arma::mat TransformTraj(arma::mat OdeTraj,arma::mat OriginLatent, List Filter_NC);


List Traj_sim_ezG2(arma::vec initial, arma::vec times, arma::vec param,
                   int gridsize,arma::vec x_r,arma::ivec x_i,double t_correct,
                   std::string transP,std::string model,
                   std::string transX);


List Traj_sim_ezG_NC(arma::vec initial, arma::vec times, arma::vec param,
                     int gridsize,arma::vec x_r,arma::ivec x_i,double t_correct,
                     std::string transP,std::string model,
                     std::string transX);


double coal_loglik3(List init, arma::mat f1, double t_correct, double lambda, int Index, std::string transX);


double volz_loglik_nh2(List init, arma::mat f1, arma::vec betaN, double t_correct, arma::ivec index,
                       std::string transX);

List ESlice_general_NC(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                       List init, arma::vec betaN, double t_correct, double lambda,
                       double coal_log, int gridsize, bool volz, std::string model,
                       std::string transX);


arma::mat ESlice_general2(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                          List init, arma::vec betaN, double t_correct, double lambda,
                          int reps, int gridsize, bool volz, std::string model,
                          std::string transX);

#endif /* LNA_FUNCTIONAL_H_*/
