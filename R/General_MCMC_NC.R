updateR0_general_NC = function(MCMC_obj, MCMC_setting, i){
  p = MCMC_obj$p
  R0_new = MCMC_obj$par[p+1] + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)

  if(R0_new <1 || R0_new > 10){
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  } # reject if out of bound

  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+p)]
  param_new[1] = R0_new

  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = "SEIR")

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                    MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = "SEIR")

  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)


  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaTs(param_new, Ode_Traj_coarse_new[,1],MCMC_setting$x_r,MCMC_setting$x_i),
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog
  }
  AR = 0
  #a = logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog
  if(is.na(a)){

    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0")
      # print(LatentTraj_new)
    }else{
      print("NA Trajectory")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+1] = R0_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}






update_ChangePoint_general_NC = function(MCMC_obj, MCMC_setting, i){
  p = MCMC_obj$p
  chpxid =  (MCMC_setting$x_i[2]+p+1):(sum(MCMC_setting$x_i[1:2])+p)
  prob = 0
  for(i in chpxid){
    newpara = MCMC_obj$par[c((p + 1):(MCMC_setting$x_i[2] + p),chpxid)]
    newpara[i-p] = MCMC_obj$par[i] * exp(runif(1,-0.1,0.1))

    Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,newpara,MCMC_setting$x_r,MCMC_setting$x_i, model = "SEIR")

    Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

    FT_new = KF_param_chol(Ode_Traj_thin_new,newpara,
                      MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = "SEIR")

    LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)
    logMultiNorm_new = log_like_traj_general2(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)


    if(MCMC_setting$likelihood == "volz"){
      coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                    betaTs(newpara,Ode_Traj_coarse_new[,1],MCMC_setting$x_r,MCMC_setting$x_i),
                                    MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
    }else{
      coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                                MCMC_obj$par[5],MCMC_setting$gridsize)
    }

    if (is.nan(logMultiNorm_new)){
      a = -Inf
    }else{
      a = dlnorm(newpara[i-p],0,MCMC_setting$chpr,log = T) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog - dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr,log = T)
      #a = dlnorm(newpara[i-2],0,MCMC_setting$chpr) + coalLog_new - MCMC_obj$coalLog - dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr)
    }

    AR = 0
    if(is.na(a)){
      AR = 0
      if(is.na(coalLog_new)){
        print("NA appears in likelihood changepoint")
        # print(LatentTraj_new)
      }else{
        print("NA appears in changpoint update trajectory")
      }
    }else if(log(runif(1,0,1)) < a) {
      AR = 1
      MCMC_obj$FT = FT_new
      MCMC_obj$par[i] = newpara[i-p]
      MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
      MCMC_obj$logMultiNorm = logMultiNorm_new
      MCMC_obj$coalLog = coalLog_new
      MCMC_obj$LatentTraj = LatentTraj_new
      prob = prob + dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr,log = T)
    }
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR,prob = prob))
}


updateTraj_general_NC = function(MCMC_obj,MCMC_setting,i){

  Res = ESlice_general_NC(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,
                                        MCMC_obj$FT, MCMC_obj$par[1:MCMC_obj$p], MCMC_setting$Init,
                                        betaN = betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                        MCMC_setting$t_correct,lambda = 1,
                                        coal_log = MCMC_obj$coalLog, MCMC_setting$gridsize,
                                        volz = (MCMC_setting$likelihood == "volz"), model = MCMC_setting$model)

  MCMC_obj$LatentTraj = Res$LatentTraj
  MCMC_obj$OriginTraj = Res$OriginTraj

  if(MCMC_setting$likelihood == "volz"){
    MCMC_obj$coalLog = volz_loglik_nh2(MCMC_setting$Init, MCMC_obj$LatentTraj,
                                       betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                       MCMC_setting$t_correct,
                                       index = MCMC_setting$x_i[3:4])
  }else{
    MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj ),MCMC_setting$t_correct,
                                   MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  return(list(MCMC_obj = MCMC_obj))
}



MCMC_setup_general2 = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,
                               a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3,c2 = 0.4,d1 = 200, d2 =40, e1 = -2.8, e2 = 0.5,chpr = 1,
                               pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18, control = list(), likelihood = "volz",model = "SEIR"){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid, t_correct)
  if(model == "SEIR"){
    p1 = 3
    p2 = 3
    x_i = c(length(changetime),3,0,2)
  }
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize = gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = x_i,
                      a1 = a1, a2 = a2,b1 = b1, b2 = b2, c1 = c1, c2 = c2,d1 = d1, d2 = d2,e1 = e1, e2 = e2,chpr = chpr,
                      pa = pa, ps1 = ps1, ps2 = ps2, pga = pga, pA = pA,control = control,
                      reps = 1, likelihood = likelihood, model = model,p=p1)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}


MCMC_initialize_general_NC = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){

  logMultiNorm = NaN
  coalLog = NaN
  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
    # print(MCMC_setting$control)
    if(is.null(MCMC_setting$control$alpha)){
      alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
      alpha2 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
    }else{
      alpha1 = MCMC_setting$control$alpha[1]
      alpha2 = MCMC_setting$control$alpha[2]
    }

    S = MCMC_setting$x_r[1] / (alpha1 + alpha2  + 1)
    E = MCMC_setting$x_r[1] * alpha1 / (alpha1 + alpha2 + 1)
    I = MCMC_setting$x_r[1] * alpha2 / (alpha1 + alpha2 + 1)

    state = c(X = S, Y = E, Z = I)

    if(is.null(MCMC_setting$control$R0)){
      R0 = runif(1,1,7)
    }else{
      R0 = MCMC_setting$control$R0
    }


    if(is.null(MCMC_setting$control$mu)){
      theta3 = exp(rnorm(1,MCMC_setting$d1,MCMC_setting$d2))
    }else{
      theta3 = MCMC_setting$control$mu
    }


    if(is.null(MCMC_setting$control$gamma)){
      gamma = exp(rnorm(1,MCMC_setting$c1, MCMC_setting$c2))
    }else{
      gamma = MCMC_setting$control$gamma
    }



    #    if(is.null(MCMC_setting$control$A)){
    #      theta4 = runif(1,0,1)
    #    }else{
    #      theta4 = MCMC_setting$control$A
    #    }

    # if(is.null(MCMC_setting$control$lambda)){
    #    lambda = rnorm(1,MCMC_setting$d1,MCMC_setting$d2)
    # }else{
    #    lambda = MCMC_setting$control$lambda
    # }

    if(is.null(MCMC_setting$control$ch)){
      ch = rlnorm(MCMC_setting$x_i[1],0,MCMC_setting$chpr)
    }else{
      ch = MCMC_setting$control$ch
    }
    if(MCMC_setting$model == "SEIR"){
      param = c(R0,theta3,gamma,ch)
    }else{
      param = c(R0, gamma, theta3, ch)
    }

    # print(param)
    #print(state)
    Ode_Traj_thin = ODE_rk45(state, MCMC_setting$times, param, MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]

    FT = KF_param_chol(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)



    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_general_noncentral(Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      #LatentTraj = Latent$SimuTraj
      LatentTraj = Latent$SimuTraj
      OriginTraj = Latent$OriginTraj
      logMultiNorm = Latent$logMultiNorm
      logOrigin = Latent$logOrigin
    }else{
      LatentTraj = MCMC_setting$control$traj
      if( sum(abs(LatentTraj[1,c(2,3)]) - c(S,I)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_traj_general_adjust(LatentTraj,Ode_Traj_coarse,
                                            FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    if(MCMC_setting$likelihood == "volz"){
      coalLog = volz_loglik_nh2(MCMC_setting$Init, LatentTraj,
                                betaTs(param, Ode_Traj_coarse[,1],MCMC_setting$x_r,MCMC_setting$x_i),
                                MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
      #print(MCMC_setting$t_correct)
    }else{
      coalLog = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj),MCMC_setting$t_correct,
                           lambda,MCMC_setting$gridsize)
    }

    print(coalLog)

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog<= - 100000000,NaN, coalLog)
    }
    plot(LatentTraj[,1],LatentTraj[,4],type="l")
  }
  LogAlpha1 = dlnorm(alpha1,MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogMu = dlnorm(theta3, MCMC_setting$e1, MCMC_setting$e2, log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dlnorm(gamma,MCMC_setting$c1,MCMC_setting$c2,log = T)
  #LogLambda = dnorm(lambda,MCMC_setting$d1,MCMC_setting$d2,log = T)

  #print(log_like_trajSIR_BD(LatentTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize,90))
  #print()
  #plot(Ode_Traj_coarse[,3])
  plot(LatentTraj[,1],LatentTraj[,3],type="l")
  paras =  c(S,E,I,param)
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,p = MCMC_setting$p,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog, OriginTraj = OriginTraj,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  print(paste("size = ", MCMC_setting$N))
  print(paste("S0 = ",S," I0 = ", I))
  print(paste("R0 = ",R0," gamma = ", gamma," mu = ", theta3))
  #print(paste("mu = ", theta3, " A = ", theta4))
  return(MCMC_obj)
}



