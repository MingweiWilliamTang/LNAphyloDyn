#include "basic_util.h"
#include "SIR_LNA.h"
#include "SIR_log_LNA.h"

#ifndef SIR_PHYLODYN_H_
#define SIR_PHYLODYN_H_
using namespace Rcpp;
using namespace arma;

double log_like_traj(arma::mat SdeTraj,arma::mat OdeTraj, List Filter,
                     int gridsize,double t_correct);

List SIR_log_KOM_Filter2(arma::mat OdeTraj,double theta1,double theta2,
                         int gridsize,std::string funname);

arma::mat ODE(arma::vec initial, arma::vec t, arma::vec param, std::string funname  = "standard");

List Traj_sim(arma::vec initial, arma::mat OdeTraj, List Filter,double t_correct);

List Traj_sim_ez(arma::vec initial, arma::vec times,double theta1, double theta2,
                 int gridsize,double t_correct = 90,std::string funname = "standard");

double coal_loglik2(List init, arma::mat f1, double t_correct, double lambda, int gridsize);

double coal_loglik(List init, arma::mat f1, double t_correct, double lambda, int gridsize);

arma::mat ESlice(arma::mat f_cur, arma::mat OdeTraj, List FTs, arma::vec state,
                 List init, double t_correct, double lambda=10, int reps=1,
                 int gridsize = 100,std::string funname  = "standard");
double volz_loglik(List init, arma::mat f1, double t_correct, double beta, int gridsize = 1);
double volz_loglik_nh(List init, arma::mat f1, arma::vec betaN, double t_correct, int gridsize = 1);
#endif /* SIR_PHYLODYN_H_ */

