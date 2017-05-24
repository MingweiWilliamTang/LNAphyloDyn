updateAlphas_log = function(MCMC_obj,MCMC_setting,i){
  #alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,- MCMC_setting$pa, MCMC_setting$pa))

  state_new = c(X = MCMC_setting$N * alpha1_new/(alpha1_new + 1),
                Y = MCMC_setting$N/(alpha1_new + 1))

  Ode_Traj_thin_new <- ODE(log(state_new),MCMC_setting$times,
                           MCMC_obj$par[3:4],"log")

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

   FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,"log")

    LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
   logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

   coalLog_new = volz_loglik(MCMC_setting$Init,LatentTraj_new, MCMC_setting$t_correct,beta =  MCMC_obj$par[3] * MCMC_setting$N,MCMC_setting$gridsize)

  if (is.nan(logMultiNorm_new)) {
    logMultiNorm_new = -Inf
    # countInf = countInf + 1
  }
  a = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T) - log(1+alpha1_new) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
    logMultiNorm_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
    ( MCMC_obj$LogAlpha1) + log(1+alpha1)

  #print(logMultiNorm_new-log_like_traj2(MCMC_obj$LatentTraj,MCMC_setting$times,log(state_new),MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,MCMC_setting$t_correct ))

  #print(logMultiNorm_new - MCMC_obj$logMultiNorm)
  if(is.na(a)){a = - Inf}
  # print(c(logMultiNorm_new,MCMC_obj$logMultiNorm,dgamma(log(alpha1_new),60,MCMC_setting$a1,log = T), dgamma(log(alpha2_new),60,MCMC_setting$a2,log = T)))
  AR = 0
  if (log(runif(1,0,1)) < a) {
    AR = 1
    state = state_new
    MCMC_obj$par[1] = state_new[1]
    MCMC_obj$par[2] = state_new[2]
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LogAlpha1 = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T)
    #    MCMC_obj$LogAlpha2 = dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



updateS1S2_log = function(MCMC_obj, MCMC_setting, i){
    s1 = MCMC_obj$par[3] / MCMC_obj$par[4] * MCMC_setting$N
   s1_new = s1 + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(s1_new <1 || s1_new > 100){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  s2_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  theta1_new =  s1_new * s2_new / MCMC_setting$N
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new)
  #   print(param_new)
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)

  Ode_Traj_thin_new <- ODE(log(MCMC_obj$par[1:2]),MCMC_setting$times,
                           param_new,"log")

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2_new,MCMC_setting$gridsize,"log")

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])

  logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  #coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  ####### volz
  coalLog_new = volz_loglik(MCMC_setting$Init,LatentTraj_new, MCMC_setting$t_correct,beta = theta1_new * MCMC_setting$N,MCMC_setting$gridsize)


  if(is.nan(logMultiNorm_new)){
    a = -Inf
  }else{
    # a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog)
    # ),1))
    a = logMultiNorm_new + log(s2_new) + dlnorm(s2_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + coalLog_new -
     MCMC_obj$logMultiNorm - MCMC_obj$coalLog - MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    #print("NA appears when update s1")
    #print(s1)
  }else if (log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = theta1_new
    MCMC_obj$par[4] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(s2_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}


updateTraj_log = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  MCMC_obj$LatentTraj = ESlice2_log(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
                                MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize,volz = T,beta = MCMC_obj$par[3] * MCMC_setting$N)
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)

  #MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  ######## volz
  MCMC_obj$coalLog = volz_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,
                                 MCMC_setting$t_correct,beta = MCMC_obj$par[3]*MCMC_setting$N,MCMC_setting$gridsize)
  ######
  # print(MCMC_obj$coalLog)
  return(list(MCMC_obj=MCMC_obj))
}


MCMC_setup_log = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5,
                      a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3, c2 = 0.2,d1 = 200, d2 =40,
                      pa = 0.1, ps1 = 0.25, ps2 = 0.5,control = list()){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid, t_correct)
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,N = N,
                      gridsize=gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,
                      a1 = a1, a2 = a2,b1 =b1, b2 = b2, c1= c1, c2 = c2,d1 = d1, d2 = d2,
                      pa = pa, ps1 = ps1, ps2 = ps2,control = control,
                      reps=1)
  cat("MCMC set up ready \n")
  return(MCMC_setting)
}


MCMC_initialize_log = function(MCMC_setting){
  #, prior_par = c(10,20,-2.3,200,40)){
  #alpha1 = exp(rgamma(1,MCMC_setting$b1,MCMC_setting$a1))
  #alpha2 = exp(rgamma(1,MCMC_setting$b2,MCMC_setting$a2))
  #S = MCMC_setting$N * alpha1 / (alpha1 + alpha2 + 1)
  #I = MCMC_setting$N * alpha2 / (alpha1 + alpha2 + 1)

  logMultiNorm = NaN
  coalLog = NaN


  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){

    #######################################
    print(MCMC_setting$control)

     if(is.null(MCMC_setting$control$alpha)){
      alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
    }else{
      alpha1 = MCMC_setting$control$alpha
    }

    S = MCMC_setting$N * alpha1 / (alpha1  + 1)
    I = MCMC_setting$N / (alpha1 + 1)
    state = c(X = S, Y = I)

    if(is.null(MCMC_setting$control$s1)){
      s1 = runif(1,1,10)
    }else{
      s1 = MCMC_setting$control$s1
    }

    if(is.null(MCMC_setting$control$s2)){
      s2 = exp(rnorm(1,MCMC_setting$c1, MCMC_setting$c2))
    }else{
      s2 = MCMC_setting$control$s2
    }

    theta2 = s2
    theta1 = s1 * s2 / MCMC_setting$N
    param = c(theta1 = theta1, theta2 = theta2)

    Ode_Traj_thin = ODE(log(state), MCMC_setting$times,
                        param,"log")
    plot(Ode_Traj_thin[,1],Ode_Traj_thin[,3],type="l")
    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]

    FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize,"log")

    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim(state,Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      LatentTraj = Latent$SimuTraj
      logMultiNorm = Latent$loglike
    }else{
      LatentTraj = MCMC_setting$control$traj
      if( sum(abs(LatentTraj[1,c(2,3)]) - c(S,I)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_traj(LatentTraj,Ode_Traj_coarse,
                                   FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    coalLog = volz_loglik(MCMC_setting$Init,LatentTraj, MCMC_setting$t_correct,
                          betaN = theta1 * MCMC_setting$N ,MCMC_setting$gridsize)

    if(coalLog <= -100000 || is.nan(coalLog)){
      coalLog = NaN
    }
    lambda = 500 # no use
  }
  print(coalLog)
  LogAlpha1 = dlnorm(alpha1,MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dlnorm(s2,MCMC_setting$c1,MCMC_setting$c2,log = T)
  LogLambda = dgamma(log(lambda),MCMC_setting$d1,MCMC_setting$d2,log = T)

  MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2, LogLambda = LogLambda)

  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  print(paste("size = ", MCMC_setting$N))
  print(paste("S0 = ",S," I0 = ", I))
  print(paste("R0 = ",s1," gamma = ", s2, " beta = ", theta1))
  return(MCMC_obj)
}



SIR_LNA_MCMC_log = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,
                                 a1 = 10, a2 = 20, b1 = 60 , b2 = 60, c1 = -2.3, c2 = 0.4, d1 = 200, d2 = 40,
                                 pa = 0.1, ps1 = 0.25, ps2 = 0.5, control = list(), updateVec = c(1,1,1,1)){
  MCMC_setting = MCMC_setup_log(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,
                            a1, a2,b1,b2,c1,c2,d1, d2,
                            pa , ps1 , ps2,control = control)
  MCMC_obj = MCMC_initialize_log(MCMC_setting)
  #print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
  # MCMC_obj$LogS2 =
  # MCMC_obj$LogLambda
  # MCMC_obj$par[1] = 8000
  # MCMC_obj$par[2] = 3000
  # MCMC_obj$par[4] = 0.05
  params = matrix(nrow = niter, ncol = 5)
  l = numeric(niter)
  l1 = l
  l2 = l
  #l3 = l
  tjs = NULL
  ARMS = matrix(nrow = niter, ncol = 4)
  for (i in 1:MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
    }
    ARvec = numeric(4)
    if(updateVec[1] == 1){
      step1 = updateAlphas_log(MCMC_obj,MCMC_setting,i)
      ARvec[1] = step1$AR
      ARvec[2] = step1$AR
      MCMC_obj = step1$MCMC_obj
    }
   # if(updateVec[2] == 1){
  #    step2 = updateS1(MCMC_obj,MCMC_setting,i)
  #    ARvec[3] = step2$AR
   #   MCMC_obj = step2$MCMC_obj
  #  }
    if(updateVec[2] == 1 && updateVec[3] == 1){
      step3 = updateS1S2_log(MCMC_obj,MCMC_setting,i)
      ARvec[4] = step3$AR
      MCMC_obj = step3$MCMC_obj
    }
    if(updateVec[4] == 1){
      MCMC_obj = updateTraj_log(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }
    tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$logMultiNorm #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = MCMC_obj$logMultiNorm + MCMC_obj$coalLog + MCMC_obj$LogAlpha1 + MCMC_obj$LogS2 -
      log(1+MCMC_obj$par[1] / MCMC_obj$par[2]) + log(MCMC_obj$par[4])

    # l3[i] = MCMC_obj$LogAlpha1 + MCMC_obj$logMultiNorm + MCMC_obj$LogS2 + MCMC_obj$coalLog + MCMC_obj$LogLambda

    ARMS[i,] = ARvec
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2,AR = ARMS))
}


