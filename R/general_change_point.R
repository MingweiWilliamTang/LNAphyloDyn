# SIRLNA_Period
updateAlphas_General = function(MCMC_obj,MCMC_setting,i){
  #alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,- MCMC_setting$pa, MCMC_setting$pa))

  #alpha2_new = alpha2 * exp(runif(x1,-0.2,0.2))
  #state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + alpha2_new + 1),
  #             Y = MCMC_setting$N * alpha2_new / (alpha1_new + alpha2_new + 1))
  state_new = c(X = MCMC_setting$x_r[1] * alpha1_new/(alpha1_new + 1),
                Y = MCMC_setting$x_r[1] / (alpha1_new + 1))

  Ode_Traj_thin_new <- General_ODE_rk45(state_new,MCMC_setting$times,MCMC_obj$par[3:(sum(MCMC_setting$x_i)+2)],
                                        MCMC_setting$x_r,MCMC_setting$x_i)


  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)],
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)
 # FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,MCMC_obj$par[3:6],MCMC_setting$gridsize,MCMC_setting$N,period = MCMC_setting$period)


  LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                          Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  #logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)
  }
  if (is.nan(logMultiNorm_new)) {
    logMultiNorm_new = -Inf
    # countInf = countInf + 1
  }
  a = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T) - 2*log(1+alpha1_new) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
    logMultiNorm_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
    ( MCMC_obj$LogAlpha1) +  2*log(1+alpha1)

  #print(logMultiNorm_new-log_like_traj2(MCMC_obj$LatentTraj,MCMC_setting$times,log(state_new),MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,MCMC_setting$t_correct ))

  #print(logMultiNorm_new - MCMC_obj$logMultiNorm)
  if(is.na(a)){a = - Inf}
  # print(c(logMultiNorm_new,MCMC_obj$logMultiNorm,dgamma(log(alpha1_new),60,MCMC_setting$a1,log = T), dgamma(log(alpha2_new),60,MCMC_setting$a2,log = T)))
  AR = 0
  if (log(runif(1,0,1)) < a) {
    AR = 1
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



##############
updateR0gamma_general = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  R0_new = MCMC_obj$par[3] + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(R0_new <1 || R0_new > 10){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  gamma_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  param_new = c(R0_new, gamma_new, MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + log(gamma_new) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
      MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0 gamma")
      # print(LatentTraj_new)
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = R0_new
    MCMC_obj$par[4] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}
############
updateR0_general = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  R0_new = MCMC_obj$par[3] + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(R0_new <1 || R0_new > 10){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  #gamma_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  param_new = c(R0_new, MCMC_obj$par[4], MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0")
      # print(LatentTraj_new)
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = R0_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}





##################
updategamma_general = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  gamma_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  param_new = c(MCMC_obj$par[3], gamma_new, MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + log(gamma_new) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
      MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
      # print(LatentTraj_new)
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[4] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



###############
updateLambdaGeneral = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[5] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = coalLog_new - MCMC_obj$coalLog + dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T) -
    MCMC_obj$LogLambda
  AR = 0

  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)
  if(is.nan(a)){
    print(coalLog_new)
  }else if(log(runif(1,0,1)) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[5] = lambda_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}


update_ChangePoint_general = function(MCMC_obj, MCMC_setting, i){
  chpxid =  (MCMC_setting$x_i[2]+3):(sum(MCMC_setting$x_i)+2)
  for( i in chpxid){
    newpara = MCMC_obj$par[c(3,4,5,chpxid)]
    newpara[i-2] = MCMC_obj$par[i] * exp(runif(1,-0.1,0.1))

    Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,newpara,MCMC_setting$x_r,MCMC_setting$x_i)
    Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

    FT_new = General_KOM_Filter(Ode_Traj_thin_new,newpara,
                                MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

    LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                             Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    if(MCMC_setting$likelihood == "volz"){
      coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                   betafs(Ode_Traj_coarse_new[,1],newpara, MCMC_setting$x_r,MCMC_setting$x_i),
                                   MCMC_setting$t_correct,
                                   MCMC_setting$gridsize)
    }else{
      coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                                MCMC_obj$par[5],MCMC_setting$gridsize)
    }

    if (is.nan(logMultiNorm_new)) {
      a = -Inf
    }else{
      a = dlnorm(newpara[i-2],0,MCMC_setting$chpr) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog + dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr)
    }
    # print(theta2_new)
    AR = 0
    if(is.na(a)){
      AR = 0
      if(is.na(coalLog_new)){
        print("NA appears in likelihood changepoint")
        # print(LatentTraj_new)
      }else{
        print("NA appears in changpoint update")
      }
    }else if(log(runif(1,0,1)) < a) {
      AR = 1
      MCMC_obj$par[i] = newpara[i-2]
      MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
      MCMC_obj$logMultiNorm = logMultiNorm_new
      MCMC_obj$FT = FT_new
      MCMC_obj$coalLog = coalLog_new
      MCMC_obj$LatentTraj = LatentTraj_new
    }
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

updateTraj_general = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  #print(MCMC_obj$par[3])
  #print(betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],MCMC_obj$LatentTraj[,1]))
  MCMC_obj$LatentTraj = ESlice_general(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,
                                      MCMC_obj$FT,MCMC_obj$par[1:2], MCMC_setting$Init,betaN = betafs(MCMC_obj$LatentTraj[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                      MCMC_setting$t_correct,lambda = MCMC_obj$par[5],
                                      reps = MCMC_setting$reps,MCMC_setting$gridsize,
                                      volz = (MCMC_setting$likelihood == "volz"))
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_traj_general(MCMC_obj$LatentTraj, MCMC_obj$Ode_Traj_coarse,
                                                MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    MCMC_obj$coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(MCMC_obj$LatentTraj),
                                      betafs(MCMC_obj$LatentTraj[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                      MCMC_setting$t_correct,
                                      MCMC_setting$gridsize)
  }else{
    MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj ),MCMC_setting$t_correct,
                                   MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  return(list(MCMC_obj=MCMC_obj))
}




MCMC_setup_general = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,
                      a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3,c2 = 0.4,d1 = 200, d2 =40, e1 = -2.8, e2 = 0.5,chpr = 1,
                      pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18, control = list(), likelihood = "volz"){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid, t_correct)
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize=gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = c(length(changetime),3),
                      a1 = a1, a2 = a2,b1 =b1, b2 = b2, c1= c1, c2 = c2,d1 = d1, d2 = d2,e1 = e1, e2 = e2,chpr = chpr,
                      pa = pa, ps1 = ps1, ps2 = ps2, pga = pga, pA = pA,control = control,
                      reps=1, likelihood = likelihood)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}


MCMC_initialize_general = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){

  logMultiNorm = NaN
  coalLog = NaN
  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
    # print(MCMC_setting$control)
    if(is.null(MCMC_setting$control$alpha)){
      alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
    }else{
      alpha1 = MCMC_setting$control$alpha
    }

    S = MCMC_setting$x_r[1] * alpha1 / (alpha1  + 1)
    I = MCMC_setting$x_r[1] / (alpha1 + 1)
    state = c(X = S, Y = I)

    if(is.null(MCMC_setting$control$R0)){
      R0 = runif(1,1,7)
    }else{
      R0 = MCMC_setting$control$R0
    }

    if(is.null(MCMC_setting$control$gamma)){
      gamma = exp(rnorm(1,MCMC_setting$c1, MCMC_setting$c2))
    }else{
      gamma = MCMC_setting$control$gamma
    }

#    if(is.null(MCMC_setting$control$mu)){
#      theta3 = exp(rnorm(1,MCMC_setting$e1,MCMC_setting$e2))
#    }else{
#      theta3 = MCMC_setting$control$mu
#    }

#    if(is.null(MCMC_setting$control$A)){
#      theta4 = runif(1,0,1)
#    }else{
#      theta4 = MCMC_setting$control$A
#    }

    if(is.null(MCMC_setting$control$lambda)){
      lambda = rnorm(1,MCMC_setting$d1,MCMC_setting$d2)
    }else{
      lambda = MCMC_setting$control$lambda
    }

    if(is.null(MCMC_setting$control$ch)){
      ch = rlnorm(MCMC_setting$x_i[1],0,MCMC_setting$chpr)
    }else{
      ch = MCMC_setting$control$ch
    }

    param = c(R0, gamma, lambda, ch)
    # print(param)
    #print(state)
    Ode_Traj_thin = General_ODE_rk45(state, MCMC_setting$times, param, MCMC_setting$x_r, MCMC_setting$x_i)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = General_KOM_Filter(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$x_r,MCMC_setting$x_i)


    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_general(Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      LatentTraj = Latent$SimuTraj
      logMultiNorm = Latent$loglike
    }else{
      LatentTraj = MCMC_setting$control$traj
      if( sum(abs(LatentTraj[1,c(2,3)]) - c(S,I)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_traj_general(LatentTraj,Ode_Traj_coarse,
                                         FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }


    if(MCMC_setting$likelihood == "volz"){
      coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj), betafs(LatentTraj[,1],param, MCMC_setting$x_r,MCMC_setting$x_i),
                               MCMC_setting$t_correct,
                               MCMC_setting$gridsize)


      #print(MCMC_setting$t_correct)
    }else{
      coalLog= coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj),MCMC_setting$t_correct,
                           lambda,MCMC_setting$gridsize)
    }

    print(coalLog)

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog<= - 100000000,NaN, coalLog)
    }
    plot(LatentTraj[,1],LatentTraj[,3],type="l")
  }
  LogAlpha1 = dlnorm(alpha1,MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogMu = dlnorm(theta3, MCMC_setting$e1, MCMC_setting$e2, log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dlnorm(gamma,MCMC_setting$c1,MCMC_setting$c2,log = T)
  LogLambda = dnorm(lambda,MCMC_setting$d1,MCMC_setting$d2,log = T)

  #print(log_like_trajSIR_BD(LatentTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize,90))
  #print()
  #plot(Ode_Traj_coarse[,3])
  plot(LatentTraj[,1],LatentTraj[,3],type="l")
  paras =  c(S,I,param)
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2,LogLambda = LogLambda)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  print(paste("size = ", MCMC_setting$N))
  print(paste("S0 = ",S," I0 = ", I))
  print(paste("R0 = ",R0," gamma = ", gamma," lambda = ", lambda))
  #print(paste("mu = ", theta3, " A = ", theta4))
  return(MCMC_obj)
}


SIR_general_MCMC = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,changetime,
                           a1 = 10, a2 = 20, b1 = 60 , b2 = 60, c1=-2.3,c2 = 0.4,d1 = 250, d2 =40, e1 = -2.8, e2 = 0.5,chpr = 1,
                           pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18, control = list(), updateVec = c(1,1,1,1,1), likelihood = "volz"){
  MCMC_setting = MCMC_setup_general(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,changetime = changetime,
                            a1, a2,b1,b2,c1,c2,d1, d2,e1,e2,chpr,
                            pa,ps1,ps2,pga,pA,control = control,likelihood = likelihood)

  MCMC_obj = MCMC_initialize_general(MCMC_setting)
  if(MCMC_setting$likelihood == "volz"){
    params = matrix(nrow = niter, ncol = 6)
    ARMS = matrix(nrow = niter, ncol = 6)
  }else{
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + 2)
    ARMS = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + 2)
  }
  l = numeric(niter)
  l1 = l
  l2 = l
  l3 = l
  tjs = NULL
  for (i in 1:MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
      print(MCMC_obj$par)
      plot(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,3],type="l")
    }
    ARvec = numeric(dim(ARMS)[2])
    if(updateVec[1] == 1){
      step1 = updateAlphas_General(MCMC_obj,MCMC_setting,i)
      # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
      MCMC_obj = step1$MCMC_obj
      ARvec[1] = step1$AR
    }
    if(updateVec[2] == 1){
      step2 = updateR0_general(MCMC_obj,MCMC_setting,i)
      # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
      MCMC_obj = step2$MCMC_obj
      ARvec[2] = step2$AR
    }
    if(updateVec[3] == 1){
      step3 = updategamma_general(MCMC_obj,MCMC_setting,i)
      # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
      MCMC_obj = step3$MCMC_obj
      ARvec[3] = step3$AR
    }

    # if(updateVec[3] == 1 && updateVec[2] == 1){
    #  stepJoint = updateR0gamma_general(MCMC_obj,MCMC_setting,i)
    #  ARvec[4] = stepJoint$AR
    #  MCMC_obj = stepJoint$MCMC_obj
    #}

    if(updateVec[4] == 1){
      steplambda = updateLambdaGeneral(MCMC_obj,MCMC_setting,i)
      ARvec[5] = steplambda$AR
      MCMC_obj = steplambda$MCMC_obj
    }

    if(updateVec[5] == 1){
      MCMC_obj = update_ChangePoint_general(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }

    if(updateVec[6] == 1){
      MCMC_obj = updateTraj_general(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }

    #step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
    #ARvec[8] = step4$AR
    #MCMC_obj = step4$MCMC_obj
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



