# SIRLNA_Period


updateS1_SIRS = function(MCMC_obj, MCMC_setting, i){
  s1 = MCMC_obj$par[4] / MCMC_obj$par[5] * MCMC_setting$N
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  s1_new = s1 + runif(1,-0.25,0.25)
  if(s1_new <1 || s1_new > 5){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  theta1_new = s1_new * MCMC_obj$par[5] / MCMC_setting$N
  param_new = c(theta1_new, MCMC_obj$par[5:7])
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize)

    LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                             Ode_Traj_coarse_new[,2:4])
    logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog)
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update s1")
    print(s1_new)
    print(MCMC_obj$par[4] / MCMC_obj$par[5] * MCMC_setting$N)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[4] = theta1_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new

    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}
###########################

updateS2_SIRS = function(MCMC_obj, MCMC_setting, i){
  s2_new = MCMC_obj$par[5] * exp(runif(1,-0.3,0.3))
  theta1_new = MCMC_obj$par[4] / MCMC_obj$par[5] * s2_new
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new, MCMC_obj$par[6:7])


  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,
                           param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)


  if (is.nan(logMultiNorm_new)) {
    a = -1
  }else{
    a = min(c(exp(dnorm(log(s2_new),MCMC_setting$c1,0.15,log = T) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
                    MCMC_obj$LogS2 ),1))
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears")
  }else if(runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[4] = theta1_new
    MCMC_obj$par[5] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dnorm(log(s2_new),MCMC_setting$c1,0.15,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new

  }
  #MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                  MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  #MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)

  #MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



#####################



updateReSus_SIRS = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta3 log(theta3) ~ N(-29.8,0.5)
  theta3_new = MCMC_obj$par[6] * exp(runif(1,-0.2,0.2))

  param_new = c(MCMC_obj$par[4:5], theta3_new, MCMC_obj$par[7])
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog + dnorm(theta3_new,-2.8,0.5,T) - dnorm(MCMC_obj$par[6],-2.8,0.5,T))
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update rsus")
    print(theta3_new)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[6] = theta3_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    # print(MCMC_obj$par)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



update_Scale_SIRS = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta4 theta ~ beta(2,2)
  A_new = MCMC_obj$par[7] + (runif(1,-0.18,0.18))
  if(A_new <=0 || A_new >= 1){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  param_new = c(MCMC_obj$par[4:6],A_new)
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog +  dbeta(A_new,2,2,log=T) - dbeta(MCMC_obj$par[7],2,2,log = T))
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update s1")
    print(s1)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[7] = A_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    # print(MCMC_obj$par)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}




##################

updateLambda_SIRS = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[8] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = min(c(exp(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
                  MCMC_obj$LogLambda), 1))
  AR = 0

  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)

  if(runif(1,0,1) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[8] = lambda_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}




updateTraj = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  MCMC_obj$LatentTraj = ESlice_SIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_obj$par[1:3],
                               MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[8],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)
  # print(MCMC_obj$coalLog)
  return(list(MCMC_obj=MCMC_obj))
}




MCMC_setup = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5,
                      a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3,d1 = 200, d2 =40){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid, t_correct)
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,N = N,
                      gridsize=gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,
                      a1 = a1, a2 = a2,b1 =b1, b2 = b2, c1= c1,d1 = d1, d2 = d2,
                      reps=1)
  #cat("MCMC set up ready \n")
  return(MCMC_setting)
}


MCMC_initialize2 = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){
  #alpha1 = exp(rgamma(1,MCMC_setting$b1,MCMC_setting$a1))
  #alpha2 = exp(rgamma(1,MCMC_setting$b2,MCMC_setting$a2))
  #S = MCMC_setting$N * alpha1 / (alpha1 + alpha2 + 1)
  #I = MCMC_setting$N * alpha2 / (alpha1 + alpha2 + 1)
  alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))

  #S = MCMC_setting$N * alpha1 / (alpha1  + 1)
  #I = MCMC_setting$N / (alpha1 + 1)
  S = 8500
  I = 500
  R = 10000
  state = c(X = S, Y = I, R = 1500)

  logMultiNorm = NaN
  coalLog = NaN
  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
    #print(coalLog)
    #s1 = runif(1,1,10)
    #s2 = exp(rnorm(1,MCMC_setting$c1,0.4))
    #theta1 = s1 * s2 / MCMC_setting$N
    #theta2 = s2
    theta1 = 0.00006
    theta2 = 0.5
    theta3 = 0.06
    theta4 = 0.0
    s2 = theta2
    s1 = theta1 * MCMC_setting$N / s2
    param = c(theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4)

    #Ode_Traj_thin <- ODE(log(state),times,param)

    Ode_Traj_thin = ODE2(state, MCMC_setting$times,
                        param)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIRS_KOM_Filter(Ode_Traj_thin,param,MCMC_setting$gridsize)


    Latent = Traj_sim_SIRS(state,Ode_Traj_coarse,FT)
    LatentTraj = Latent$SimuTraj
    lambda = 400
    coalLog = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj),MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
    # lambda = exp(rgamma(1,16,4))

    logMultiNorm = Latent$loglike
  }
  LogAlpha1 = dnorm(log(alpha1),MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4,log = T)
  LogLambda = dgamma(log(lambda),MCMC_setting$d1,MCMC_setting$d2,log = T)
  print(log_like_trajSIRS(LatentTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize,90))
  print(logMultiNorm)
  #plot(Ode_Traj_coarse[,3])
  MCMC_obj = list(par = c(S,I,R,theta1,theta2,theta3,theta4,lambda),LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2, LogLambda = LogLambda)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  return(MCMC_obj)
}



SIRS_LNA_MCMC_standard = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,
                                 a1 = 10, a2 = 20, b1 = 60 , b2 = 60, c1 = -2.3, d1 = 200, d2 = 40){
  MCMC_setting = MCMC_setup(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,
                            a1, a2,b1,b2,c1,d1, d2)
  MCMC_obj = MCMC_initialize2(MCMC_setting)
  #print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
  # MCMC_obj$LogS2 =
  # MCMC_obj$LogLambda
  MCMC_obj$par[1] = 8500
  MCMC_obj$par[2] = 500
  MCMC_obj$par[3] = 10000
  MCMC_obj$par[5] = 0.5
  MCMC_obj$par[6] = 0.06
  MCMC_obj$par[7] = 0
  MCMC_obj$par[8] = 400
  params = matrix(nrow = niter, ncol = 8)
  ARMS = matrix(nrow = niter, ncol = 8)
  l = numeric(niter)
  l1 = l
  l2 = l
  l3 = l
  tjs = NULL
  for (i in 1:MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
    }
    ARvec = numeric(8)
    #step1 = updateAlphas(MCMC_obj,MCMC_setting,i)
    #  print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
    #MCMC_obj = step1$MCMC_obj
    step2 = updateS1_SIRS(MCMC_obj,MCMC_setting,i)
    ARvec[4] = step2$AR
    MCMC_obj = step2$MCMC_obj
   #    stepRecover = updateS2_SIRS(MCMC_obj,MCMC_setting,i)
    #  ARvec[5] = stepRecover$AR
     #  MCMC_obj = stepRecover$MCMC_obj
      stepReS = updateReSus_SIRS(MCMC_obj,MCMC_setting,i)
      ARvec[6] = stepReS$AR
      MCMC_obj = stepReS$MCMC_obj
    #stepA = update_Scale_SIRS(MCMC_obj, MCMC_setting,i)
    #  ARvec[7] = stepA$AR
    #  MCMC_obj = stepA$MCMC_obj
    MCMC_obj = updateTraj(MCMC_obj,MCMC_setting,i)$MCMC_obj
   step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
   ARvec[8] = step4$AR
   MCMC_obj = step4$MCMC_obj
   ARMS[i,] = ARvec
    tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$logMultiNorm #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = MCMC_obj$logMultiNorm + MCMC_obj$coalLog
    l3[i] = MCMC_obj$LogAlpha1 + MCMC_obj$logMultiNorm + MCMC_obj$LogS2 + MCMC_obj$coalLog + MCMC_obj$LogLambda
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2, l3 =l3,AR = ARMS))
}


effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1 /(2 * Traj[,3] * beta / Traj[,2]))
  }else{
    return(Traj[,3] / lambda)
  }
}
