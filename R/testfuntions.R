Geweke_S1_S2_Traj = function(gridsize,nn=100){
  times = seq(0,100,by=0.01)
  a = numeric(nn)
  b = numeric(nn)
    #Traj2 = Traj_sim_ez(c(10000,1000),seq(0,100,by=0.01),0.000015,0.05,gridsize = 100,90,"standard")$Simu
    gridset = seq(1,length(times),by=gridsize)
    #coal_obj =  list(samp_times = c(seq(0,70,length.out=200)),n_sampled = c(rep(2,200)),
    #                  coal_times = G$coal_times )
    S = 10000
    I = 1000
    N = 11000
    state = c(X = S, Y = I)
    s1 = runif(1,3,6)
    s2 = 0.05
    theta1 = s1 * s2 / N
    theta2 = s2
    param = c(theta1 = theta1, theta2 = theta2)
    Tjsim = Traj_sim_ez(state,times,theta1,theta2,gridsize,t_correct = 90,"standard")
    Traj2 = Tjsim$SimuTraj
    logMultiNorm = Tjsim$loglike
    G = coalsim_thin_sir(Traj2[,c(1,3)],t_correct = 90,
                         samp_times = c(seq(0,50,length.out=200)),c(rep(3,200)),500)
    coal_obj = list(samp_times = c(seq(0,50,length.out=200)),n_sampled = c(rep(3,200)),
                    coal_times = G$coal_times )
    MCMC_setting = MCMC_setup(coal_obj, times, 90,N, gridsize,niter=1000,burn = 0,thin=1,a1=10,a2=20,b1=60,b2=60,c1=-2.3,
                              d1=200,d2=40)

    Ode_Traj_thin = ODE(state, MCMC_setting$times,
                        param,"standard")

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize,"standard")
    lambda = 500
    coalLog = coal_loglik(MCMC_setting$Init,LogTraj(Traj2),MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
    LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4,log = T)
    MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = Traj2, logMultiNorm = logMultiNorm,
                    Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                    LogAlpha1 = NULL, LogS2 = LogS2, LogLambda = NULL)
    for(i in 1:nn){
      MCMC_obj = updateS1(MCMC_obj,MCMC_setting,1)$MCMC_obj
      MCMC_obj = updateS2(MCMC_obj,MCMC_setting,1)$MCMC_obj
      if(i %% 50 == 0){
        print(paste("finish ", i, " iterations"))
      }
      #MCMC_obj = updateS2(MCMC_obj,MCMC_setting,1)$MCMC_obj
     for(j in 1){
       MCMC_obj = updateTraj(MCMC_obj, MCMC_setting,1)$MCMC_obj
     }
    # print(MCMC_obj$par)
      G = coalsim_thin_sir(MCMC_obj$LatentTraj[,c(1,3)],t_correct = 90,
                           samp_times = c(seq(0,50,length.out=200)),c(rep(3,200)),500)
      coal_obj = list(samp_times = c(seq(0,50,length.out=200)),n_sampled = c(rep(3,200)),
                      coal_times = G$coal_times )
      MCMC_setting = MCMC_setup(coal_obj, times, 90,N,gridsize,niter = 1000, burn = 0,thin=1,a1=10,a2=20,b1=60,b2=60,c1=-2.3,
                                d1=200,d2=40)
      a[i] = MCMC_obj$par[3] * 11000/ MCMC_obj$par[4]
      b[i] = MCMC_obj$par[4]
    }
    return(list(theta1 = a, theta2 = b))
}



Geweke_S1_Traj_lambda = function(gridsize,nn=100){
  times = seq(0,100,by=0.01)
  a = numeric(nn)
  b = numeric(nn)
  #Traj2 = Traj_sim_ez(c(10000,1000),seq(0,100,by=0.01),0.000015,0.05,gridsize = 100,90,"standard")$Simu
  gridset = seq(1,length(times),by=gridsize)
  #coal_obj =  list(samp_times = c(seq(0,70,length.out=200)),n_sampled = c(rep(2,200)),
  #                  coal_times = G$coal_times )
  S = 10000
  I = 1000
  N = 11000
  state = c(X = S, Y = I)
  s1 = runif(1,3,6)
  s2 = 0.05
  lambda = runif(1,300,900)
  theta1 = s1 * s2 / N
  theta2 = s2
  param = c(theta1 = theta1, theta2 = theta2)
  Tjsim = Traj_sim_ez(state,times,theta1,theta2,gridsize,t_correct = 90,"standard")
  Traj2 = Tjsim$SimuTraj
  logMultiNorm = Tjsim$loglike
  G = coalsim_thin_sir(Traj2[,c(1,3)],t_correct = 90,
                       samp_times = c(seq(0,50,length.out=200)),c(rep(3,200)),500)
  coal_obj = list(samp_times = c(seq(0,50,length.out=200)),n_sampled = c(rep(3,200)),
                  coal_times = G$coal_times )
  MCMC_setting = MCMC_setup(coal_obj, times, 90,N, gridsize,niter=1000,burn = 0,thin=1,a1=10,a2=20,b1=60,b2=60,c1=-2.3,
                            d1=200,d2=40)

  Ode_Traj_thin = ODE(state, MCMC_setting$times,
                      param,"standard")

  Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


  FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize,"standard")
  lambda = 500
  coalLog = coal_loglik(MCMC_setting$Init,LogTraj(Traj2),MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4,log = T)
  MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = Traj2, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = NULL, LogS2 = LogS2, LogLambda = NULL)
  for(i in 1:nn){
    MCMC_obj = updateS1(MCMC_obj,MCMC_setting,1)$MCMC_obj
    MCMC_obj = updateLambdaUnifProp(MCMC_obj, MCMC_setting, 1)$MCMC_obj
    if(i %% 50 == 0){
      print(paste("finish ", i, " iterations"))
    }
    #MCMC_obj = updateS2(MCMC_obj,MCMC_setting,1)$MCMC_obj
    for(j in 1){
      MCMC_obj = updateTraj(MCMC_obj, MCMC_setting,1)$MCMC_obj
    }
    # print(MCMC_obj$par)
    G = coalsim_thin_sir(MCMC_obj$LatentTraj[,c(1,3)],t_correct = 90,
                         samp_times = c(seq(0,50,length.out=200)),c(rep(3,200)),500)
    coal_obj = list(samp_times = c(seq(0,50,length.out=200)),n_sampled = c(rep(3,200)),
                    coal_times = G$coal_times )
    MCMC_setting = MCMC_setup(coal_obj, times, 90,N,gridsize,niter = 1000, burn = 0,thin=1,a1=10,a2=20,b1=60,b2=60,c1=-2.3,
                              d1=200,d2=40)
    a[i] = MCMC_obj$par[3] * 11000/ MCMC_obj$par[4]
    b[i] = MCMC_obj$par[5]
  }
  return(list(s1=a,lambda=b))
}

